#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <unordered_map>
#include <functional> // for std::hash
#include <chrono>     

// Constants
const double k = 8.99e9; 
const double q = 1.6e-19; 

struct Particle {
    double x;
    double y;
    char charge; // 'p' and 'e' 
};

// Custom hash function for std::pair<int, int>
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

// Euclidean distance
double calculateDistance(const Particle& p1, const Particle& p2) {
    return std::sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
}

// calculate the force magnitude
double calculateForceMagnitude(double distance) {
    distance *= 1e-10;
    return (k * q * q) / (distance * distance);
}

// Calculate the cell indices based on the cutoff radius
std::pair<int, int> getCellIndex(const Particle& particle, double cellSize) {
    int cellX = static_cast<int>(particle.x / cellSize);
    int cellY = static_cast<int>(particle.y / cellSize);
    return {cellX, cellY};
}

// calculate the error between computed and oracle values
double calculatePercentageError(const std::vector<double>& computed, const std::vector<double>& oracle) {
    double totalRelativeError = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < computed.size(); ++i) {
        double oracleForce = oracle[i];
        double computedForce = computed[i];

        if (std::isinf(oracleForce)) {
            continue;
        }

        
        double relativeError = std::abs(computedForce - oracleForce) / std::abs(oracleForce);
        totalRelativeError += relativeError;
        count++;
    }

    return (totalRelativeError / count) * 100.0; 
}



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <cutoff_radius>" << std::endl;
        return 1;
    }

    std::string inputFile = "dataset/particles.csv";
    std::string outputFile = "output/mode1output.csv";
    std::string oracleFile = "dataset/oracle.csv";  
    double cutoffRadius = std::stod(argv[1]);
    double cellSize = cutoffRadius; 

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

        // Parse x, y coordinates and charge
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

    std::unordered_map<std::pair<int, int>, std::vector<Particle>, pair_hash> grid;

    // Populate the grid with particles
    for (const auto& particle : particles) {
        auto cellIndex = getCellIndex(particle, cellSize);
        grid[cellIndex].push_back(particle);
    }

    std::vector<double> netForces(particles.size(), 0.0);

    
    auto calculationStart = std::chrono::high_resolution_clock::now();

    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p1 = particles[i];
        auto cellIndex = getCellIndex(p1, cellSize);

        // Check this cell and neighboring cells
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                auto neighborCell = std::make_pair(cellIndex.first + dx, cellIndex.second + dy);

                
                if (grid.count(neighborCell) > 0) {
                    for (const auto& p2 : grid[neighborCell]) {
                        if (p1.x != p2.x || p1.y != p2.y) { // Avoid self-interaction
                            double distance = calculateDistance(p1, p2);
                            if (distance <= cutoffRadius) {
                                double forceMagnitude = calculateForceMagnitude(distance);

                                
                                if (p1.charge == p2.charge) {
                                    netForces[i] += forceMagnitude; // Repulsive force
                                } else {
                                    netForces[i] -= forceMagnitude; // Attractive force
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    
    auto calculationEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> calculationTime = calculationEnd - calculationStart;
    std::cout << "Force Calculation Stage Time: " << calculationTime.count() << " seconds" << std::endl;

    double combinedTime = parsingTime.count() + calculationTime.count();
    std::cout << "Total Time: " << combinedTime << " seconds" << std::endl;

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

    std::cout << "Sequential force calculation complete with spatial partitioning. Output saved to " << outputFile << std::endl;

    return 0;
}
