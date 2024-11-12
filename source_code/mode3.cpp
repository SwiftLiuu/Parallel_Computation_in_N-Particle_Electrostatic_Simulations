#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <unordered_map>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <mpi.h>
#include <chrono>

// Constants
const double k = 8.99e9; 
const double q = 1.6e-19; 

struct Particle {
    double x;
    double y;
    char charge;
    int index;
};

// Custom hash function for grid cell coordinates
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

// Task structure for each cell
struct CellTask {
    std::pair<int, int> cellIndex;
    std::vector<Particle> particles;
};

// Global variables for thread-safe queue and condition
std::queue<CellTask> cellQueue;
std::mutex queueMutex;
std::condition_variable queueCV;
bool workDone = false;

double calculateDistance(const Particle& p1, const Particle& p2) {
    return std::sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
}

double calculateForceMagnitude(double distance) {
    distance *= 1e-10;
    return (k * q * q) / (distance * distance);
}

std::pair<int, int> getCellIndex(const Particle& particle, double cellSize) {
    int cellX = static_cast<int>(particle.x / cellSize);
    int cellY = static_cast<int>(particle.y / cellSize);
    return {cellX, cellY};
}









void workerFunction(const std::unordered_map<std::pair<int, int>, std::vector<Particle>, pair_hash>& grid,
                    std::vector<double>& netForces, double cutoffRadius) {
    while (true) {
        CellTask task;

        {
            std::unique_lock<std::mutex> lock(queueMutex);
            queueCV.wait(lock, [] { return !cellQueue.empty() || workDone; });

            if (workDone && cellQueue.empty()) break;

            task = cellQueue.front();
            cellQueue.pop();
        }

        const auto& particles = task.particles;
        const auto& cellIndex = task.cellIndex;

        for (const auto& p1 : particles) {
            double totalForce = 0.0;

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    auto neighborIndex = std::make_pair(cellIndex.first + dx, cellIndex.second + dy);

                    if (grid.count(neighborIndex) == 0) {
                        // Neighbor cell is missing, so continue to the next neighbor.
                        continue;
                    }

                    const auto& neighborParticles = grid.at(neighborIndex);
                    for (const auto& p2 : neighborParticles) {
                        if (p1.index != p2.index) {
                            double distance = calculateDistance(p1, p2);
                            if (distance <= cutoffRadius) {
                                double forceMagnitude = calculateForceMagnitude(distance);
                                totalForce += (p1.charge == p2.charge ? forceMagnitude : -forceMagnitude);
                            }
                        }
                    }
                }
            }
            netForces[p1.index] = totalForce;
        }
    }
}










void createCellTasks(const std::unordered_map<std::pair<int, int>, std::vector<Particle>, pair_hash>& grid) {
    for (const auto& cell : grid) {
        CellTask task = {cell.first, cell.second};
        std::lock_guard<std::mutex> lock(queueMutex);
        cellQueue.push(task);
    }
    workDone = true;  // Set workDone to true only after all tasks are added
    queueCV.notify_all();  // Notify all workers that tasks are ready
}






double calculatePercentageError(const std::vector<double>& computed, const std::vector<double>& oracle) {
    double totalRelativeError = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < computed.size(); ++i) {
        double oracleForce = oracle[i];
        double computedForce = computed[i];

        if (std::isinf(oracleForce)) continue;

        double relativeError = std::abs(computedForce - oracleForce) / std::abs(oracleForce);
        totalRelativeError += relativeError;
        count++;
    }
    return (totalRelativeError / count) * 100.0;
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, numLeaders;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numLeaders);

    if (argc < 3) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <cutoff_radius> <num_workers>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    double cutoffRadius = std::stod(argv[1]);
    int numWorkers = std::stoi(argv[2]);
    double cellSize = cutoffRadius;

    std::string inputFile = "dataset/particles.csv";
    std::string outputFile = "output/mode3output.csv";
    std::string oracleFile = "dataset/oracle.csv";

    std::vector<Particle> particles;
    int numParticles = 0;

    // Parsing Data Stage (only for rank 0)
    if (rank == 0) {
        auto parsingStart = std::chrono::high_resolution_clock::now();

        std::cout << "Rank " << rank << ": Reading particles from file." << std::endl;
        std::ifstream input(inputFile);
        if (!input.is_open()) {
            std::cerr << "Error opening input file: " << inputFile << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int index = 0;
        std::string line;
        while (std::getline(input, line)) {
            std::stringstream ss(line);
            std::string x_str, y_str;
            char charge;
            std::getline(ss, x_str, ',');
            std::getline(ss, y_str, ',');
            ss >> charge;

            double x = std::stod(x_str);
            double y = std::stod(y_str);
            particles.push_back({x, y, charge, index++});
        }
        input.close();
        numParticles = particles.size();

        auto parsingEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parsingTime = parsingEnd - parsingStart;
        std::cout << "Parsing Data Stage Time: " << parsingTime.count() << " seconds" << std::endl;

        std::cout << "Rank " << rank << ": Finished reading particles. Total particles: " << numParticles << std::endl;
    }

    // Broadcast the number of particles to all processes
    MPI_Bcast(&numParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize particles vector for non-root processes
    if (rank != 0) {
        particles.resize(numParticles);
    }

    // Broadcast all particles to each process
    MPI_Bcast(particles.data(), numParticles * sizeof(Particle), MPI_BYTE, 0, MPI_COMM_WORLD);

    int particlesPerLeader = numParticles / numLeaders;
    int remainder = numParticles % numLeaders;
    int start = rank * particlesPerLeader + std::min(rank, remainder);
    int end = start + particlesPerLeader + (rank < remainder ? 1 : 0);

    std::cout << "Rank " << rank << ": Handling particles " << start << " to " << end - 1 << std::endl;

    std::vector<double> netForces(numParticles, 0.0);


    auto partitionStart = std::chrono::high_resolution_clock::now();

    std::unordered_map<std::pair<int, int>, std::vector<Particle>, pair_hash> grid;
    for (const auto& particle : particles) {
        auto cellIndex = getCellIndex(particle, cellSize);
        grid[cellIndex].push_back(particle);
    }

    auto partitionEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> partitionTime = partitionEnd - partitionStart;
    std::cout << "Rank " << rank << ": Data Partition Stage Time: " << partitionTime.count() << " seconds" << std::endl;


    auto calculationStart = std::chrono::high_resolution_clock::now();

    createCellTasks(grid);

    std::vector<std::thread> workers;
    for (int i = 0; i < numWorkers; ++i) {
        workers.emplace_back(workerFunction, std::ref(grid), std::ref(netForces), cutoffRadius);
    }

    for (auto& worker : workers) {
        worker.join();
    }

    auto calculationEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> calculationTime = calculationEnd - calculationStart;
    std::cout << "Rank " << rank << ": Force Calculation Stage Time: " << calculationTime.count() << " seconds" << std::endl;


    std::vector<double> globalNetForces(numParticles);
    MPI_Gather(
        netForces.data() + start, end - start, MPI_DOUBLE,
        globalNetForces.data(), end - start, MPI_DOUBLE, 0, MPI_COMM_WORLD
    );

    if (rank == 0) {
        std::ifstream oracleInput(oracleFile);
        std::vector<double> oracleForces(numParticles);
        double force;
        int i = 0;
        while (oracleInput >> force) {
            oracleForces[i++] = force;
        }
        oracleInput.close();

        double errorPercentage = calculatePercentageError(globalNetForces, oracleForces);
        std::cout << "Error Percentage: " << errorPercentage << "%" << std::endl;

        std::ofstream output(outputFile);
        for (const auto& netForce : globalNetForces) {
            output << netForce << std::endl;
        }
        output.close();

        std::cout << "Parallel force calculation complete with " << numLeaders
                  << " leaders, each with " << numWorkers << " workers. Output saved to "
                  << outputFile << std::endl;
    }

    MPI_Finalize();
    return 0;
}
