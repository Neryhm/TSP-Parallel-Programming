#include <iostream>
#include <ctime>

template <int N>
class TSPBnB {
    static const int INF = 1000000000; // Infinity value
    static const unsigned long long LOG_INTERVAL = 1000000ULL; // Log every 1,000,000 nodes
    int costMatrix[N][N];              // Cost matrix (to be provided externally)
    int finalCost;                     // Minimum cost of the tour
    int finalPath[N + 1];              // Optimal path
    bool visited[N];                   // Tracks visited cities
    double timeTaken;                  // Total computation time
    unsigned long long nodeCount;      // Number of nodes explored
    clock_t startTime;                 // Start time for timing (declared here)

    // Copies current path to final path
    void copyToFinal(int currPath[]) {
        for (int i = 0; i < N; i++) {
            finalPath[i] = currPath[i];
        }
        finalPath[N] = currPath[0]; // Return to start
    }

    // Computes a simple lower bound for pruning
    int lowerBound(int currPath[], int level, int currCost) {
        int reducedMatrix[N][N];
        int lb = 0;
    
        // Copy current matrix and set visited cities' rows/columns to INF
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                reducedMatrix[i][j] = costMatrix[i][j];
                if (visited[i] || (j == currPath[level - 1] && level > 1)) {
                    reducedMatrix[i][j] = INF;
                }
            }
        }
    
        // Row reduction
        for (int i = 0; i < N; i++) {
            if (!visited[i]) {
                int min = INF;
                for (int j = 0; j < N; j++) {
                    if (reducedMatrix[i][j] < min) min = reducedMatrix[i][j];
                }
                if (min != INF && min > 0) {
                    lb += min;
                    for (int j = 0; j < N; j++) {
                        if (reducedMatrix[i][j] != INF) reducedMatrix[i][j] -= min;
                    }
                }
            }
        }
    
        // Column reduction
        for (int j = 0; j < N; j++) {
            int min = INF;
            for (int i = 0; i < N; i++) {
                if (reducedMatrix[i][j] < min) min = reducedMatrix[i][j];
            }
            if (min != INF && min > 0) {
                lb += min;
                for (int i = 0; i < N; i++) {
                    if (reducedMatrix[i][j] != INF) reducedMatrix[i][j] -= min;
                }
            }
        }
    
        return lb + currCost; // Add current path cost
    }

    // Recursive BnB function with progress debugging
    void tspRec(int currPath[], int level, int currCost) {
        nodeCount++;
        if (nodeCount % LOG_INTERVAL == 0) {
            double elapsed = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
            std::cout << "Nodes explored: " << nodeCount << ", Current level: " << level
                      << ", Best cost: " << finalCost << ", Time elapsed: " << elapsed << " seconds\n";
        }

        if (level == N) { // Base case: all cities visited
            int returnCost = costMatrix[currPath[level - 1]][currPath[0]];
            if (returnCost != INF) {
                int totalCost = currCost + returnCost;
                if (totalCost < finalCost) {
                    copyToFinal(currPath);
                    finalCost = totalCost;
                    std::cout << "New best solution found: " << finalCost << " at node " << nodeCount << "\n";
                }
            }
            return;
        }

        // Explore next cities
        for (int i = 0; i < N; i++) {
            if (!visited[i] && costMatrix[currPath[level - 1]][i] != INF) {
                int newCost = currCost + costMatrix[currPath[level - 1]][i];
                int lb = lowerBound(currPath, level, newCost);
                if (newCost + lb < finalCost) { // Prune if cost exceeds best
                    currPath[level] = i;
                    visited[i] = true;
                    tspRec(currPath, level + 1, newCost);
                    visited[i] = false; // Backtrack
                }
            }
        }
    }

public:
    // Constructor: Initializes with provided matrix
    TSPBnB(int (&matrix)[N][N]) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                costMatrix[i][j] = matrix[i][j];
            }
            visited[i] = false;
        }
        finalCost = INF;
    }

    // Solves the TSP
    void solve() {
        nodeCount = 0;
        startTime = clock(); // Initialize startTime here
        int currPath[N + 1];
        currPath[0] = 0; // Start at city 0
        visited[0] = true;
        tspRec(currPath, 1, 0);
        clock_t end = clock();
        timeTaken = static_cast<double>(end - startTime) / CLOCKS_PER_SEC;
    }

    // Prints the solution
    void printSolution() {
        std::cout << "Minimum cost: " << finalCost << "\nPath: ";
        for (int i = 0; i <= N; i++) {
            std::cout << finalPath[i] + 1 << " "; // 1-based indexing
        }
        std::cout << "\nTime taken: " << timeTaken << " seconds\n";
    }
};

int main() {
    const int N = 13; // Number of cities
    int matrix[N][N] = {
        {0, 45, 72, 19, 84, 36, 91, 47, 63, 28, 95, 14, 67},
        {45, 0, 13, 58, 29, 94, 15, 60, 37, 73, 21, 87, 42},
        {72, 13, 0, 81, 93, 17, 80, 34, 89, 41, 66, 22, 97},
        {19, 58, 81, 0, 45, 70, 26, 91, 38, 68, 34, 99, 52},
        {84, 29, 93, 45, 0, 61, 16, 72, 83, 49, 94, 20, 65},
        {36, 94, 17, 70, 61, 0, 85, 41, 96, 22, 67, 33, 88},
        {91, 15, 80, 26, 16, 85, 0, 60, 45, 90, 36, 81, 47},
        {47, 60, 34, 91, 72, 41, 60, 0, 66, 12, 77, 43, 88},
        {63, 37, 89, 38, 83, 96, 45, 66, 0, 71, 27, 82, 48},
        {28, 73, 41, 68, 49, 22, 90, 12, 71, 0, 56, 21, 76},
        {95, 21, 66, 34, 94, 67, 36, 77, 27, 56, 0, 61, 17},
        {14, 87, 22, 99, 20, 33, 81, 43, 82, 21, 61, 0, 66},
        {67, 42, 97, 52, 65, 88, 47, 88, 48, 76, 17, 66, 0}
    };

    TSPBnB<N> tsp(matrix);
    tsp.solve();
    tsp.printSolution();
    return 0;
}