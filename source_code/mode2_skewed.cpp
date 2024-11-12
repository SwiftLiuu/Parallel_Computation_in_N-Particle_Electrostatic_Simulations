#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <unordered_map>
#include <pthread.h>
#include <chrono>

// Constants
const double k = 8.99e9; 
const double q = 1.6e-19; 

struct Particle {
    double x;
    double y;
    char charge; // 'p' and 'e' 
};

// Custom hash function for grid cell coordinates
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

// Thread argument structure
struct ThreadArgs {
    int start_index;
    int end_index;
    double cutoff_radius;
    const std::vector<Particle>* particles;
    std::unordered_map<std::pair<int, int>, std::vector<Particle>, pair_hash>* grid;
    std::vector<double>* netForces;
    double cellSize;
};


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

void* computeForces(void* args) {
    ThreadArgs* data = (ThreadArgs*)args;
    int start = data->start_index;
    int end = data->end_index;
    double cutoff_radius = data->cutoff_radius;
    double cellSize = data->cellSize;
    const std::vector<Particle>& particles = *data->particles;
    std::unordered_map<std::pair<int, int>, std::vector<Particle>, pair_hash>& grid = *data->grid;
    std::vector<double>& netForces = *data->netForces;

    for (int i = start; i < end; ++i) {
        const Particle& p1 = particles[i];
        auto cellIndex = getCellIndex(p1, cellSize);
        double totalForce = 0.0;

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                auto neighborCell = std::make_pair(cellIndex.first + dx, cellIndex.second + dy);

                if (grid.count(neighborCell) > 0) {
                    for (const auto& p2 : grid[neighborCell]) {
                        if (p1.x != p2.x || p1.y != p2.y) { // Avoid self-interaction
                            double distance = calculateDistance(p1, p2);
                            if (distance <= cutoff_radius) {
                                double forceMagnitude = calculateForceMagnitude(distance);
                                totalForce += (p1.charge == p2.charge ? forceMagnitude : -forceMagnitude);
                            }
                        }
                    }
                }
            }
        }
        netForces[i] = totalForce;
    }
    return nullptr;
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
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <cutoff_radius> <num_threads>" << std::endl;
        return 1;
    }

    double cutoffRadius = std::stod(argv[1]);
    int numThreads = std::stoi(argv[2]);
    double cellSize = cutoffRadius;

    std::string inputFile = "dataset/particles_skewed.csv";
    std::string outputFile = "output/output_skewed.csv";
    std::string oracleFile = "dataset/oracle_skewed.csv";

    std::vector<Particle> particles;

    
    auto parsingStart = std::chrono::high_resolution_clock::now();

    std::ifstream input(inputFile);
    if (!input.is_open()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return 1;
    }

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
        particles.push_back({x, y, charge});
    }
    input.close();

    auto parsingEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parsingTime = parsingEnd - parsingStart;
    std::cout << "Parsing Data Stage Time: " << parsingTime.count() << " seconds" << std::endl;

    int numParticles = particles.size();
    std::vector<double> netForces(numParticles, 0.0);

    
    auto partitionStart = std::chrono::high_resolution_clock::now();

    std::unordered_map<std::pair<int, int>, std::vector<Particle>, pair_hash> grid;
    for (const auto& particle : particles) {
        auto cellIndex = getCellIndex(particle, cellSize);
        grid[cellIndex].push_back(particle);
    }

    auto partitionEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> partitionTime = partitionEnd - partitionStart;
    std::cout << "Data Partition Stage Time: " << partitionTime.count() << " seconds" << std::endl;

    pthread_t threads[numThreads];
    ThreadArgs threadArgs[numThreads];
    int chunkSize = numParticles / numThreads;

    auto calculationStart = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numThreads; ++i) {
        threadArgs[i].start_index = i * chunkSize;
        threadArgs[i].end_index = (i == numThreads - 1) ? numParticles : (i + 1) * chunkSize;
        threadArgs[i].cutoff_radius = cutoffRadius;
        threadArgs[i].particles = &particles;
        threadArgs[i].grid = &grid;
        threadArgs[i].netForces = &netForces;
        threadArgs[i].cellSize = cellSize;

        if (pthread_create(&threads[i], nullptr, computeForces, &threadArgs[i]) != 0) {
            std::cerr << "Error creating thread " << i << std::endl;
            return 1;
        }
    }

    for (int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], nullptr);
    }

    auto calculationEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> calculationTime = calculationEnd - calculationStart;
    std::cout << "Force Calculation Stage Time: " << calculationTime.count() << " seconds" << std::endl;

    double combinedTime = parsingTime.count() + partitionTime.count() + calculationTime.count();
    std::cout << "Total Time: " << combinedTime << " seconds" << std:: endl;

    std::ifstream oracleInput(oracleFile);
    if (!oracleInput.is_open()) {
        std::cerr << "Error opening oracle file: " << oracleFile << std::endl;
        return 1;
    }

    std::vector<double> oracleForces;
    double force;
    while (oracleInput >> force) {
        oracleForces.push_back(force);
    }
    oracleInput.close();

    double errorPercentage = calculatePercentageError(netForces, oracleForces);
    std::cout << "Error Percentage: " << errorPercentage << "%" << std::endl;

    std::ofstream output(outputFile);
    if (!output.is_open()) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return 1;
    }

    for (const auto& netForce : netForces) {
        output << netForce << std::endl;
    }
    output.close();

    std::cout << "Parallel force calculation complete with " << numThreads << " threads. Output saved to " << outputFile << std::endl;

    return 0;
}
