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
    const int N = 12; // Number of cities


    // 12x12 Matrix
    int matrix[N][N] = {
        {0, 23, 56, 78, 12, 45, 89, 34, 67, 91, 28, 50},
        {23, 0, 71, 19, 83, 62, 47, 95, 31, 64, 76, 42},
        {56, 71, 0, 88, 25, 93, 16, 70, 52, 39, 81, 27},
        {78, 19, 88, 0, 44, 66, 92, 38, 75, 11, 63, 97},
        {12, 83, 25, 44, 0, 59, 21, 86, 48, 72, 35, 90},
        {45, 62, 93, 66, 59, 0, 74, 29, 82, 17, 53, 41},
        {89, 47, 16, 92, 21, 74, 0, 65, 33, 87, 26, 68},
        {34, 95, 70, 38, 86, 29, 65, 0, 51, 94, 22, 77},
        {67, 31, 52, 75, 48, 82, 33, 51, 0, 69, 84, 15},
        {91, 64, 39, 11, 72, 17, 87, 94, 69, 0, 46, 58},
        {28, 76, 81, 63, 35, 53, 26, 22, 84, 46, 0, 73},
        {50, 42, 27, 97, 90, 41, 68, 77, 15, 58, 73, 0}
    };


    // // 13x13 Matrix
    // int matrix[N][N] = {
    //     {0, 45, 72, 19, 84, 36, 91, 47, 63, 28, 95, 14, 67},
    //     {45, 0, 13, 58, 29, 94, 15, 60, 37, 73, 21, 87, 42},
    //     {72, 13, 0, 81, 93, 17, 80, 34, 89, 41, 66, 22, 97},
    //     {19, 58, 81, 0, 45, 70, 26, 91, 38, 68, 34, 99, 52},
    //     {84, 29, 93, 45, 0, 61, 16, 72, 83, 49, 94, 20, 65},
    //     {36, 94, 17, 70, 61, 0, 85, 41, 96, 22, 67, 33, 88},
    //     {91, 15, 80, 26, 16, 85, 0, 60, 45, 90, 36, 81, 47},
    //     {47, 60, 34, 91, 72, 41, 60, 0, 66, 12, 77, 43, 88},
    //     {63, 37, 89, 38, 83, 96, 45, 66, 0, 71, 27, 82, 48},
    //     {28, 73, 41, 68, 49, 22, 90, 12, 71, 0, 56, 21, 76},
    //     {95, 21, 66, 34, 94, 67, 36, 77, 27, 56, 0, 61, 17},
    //     {14, 87, 22, 99, 20, 33, 81, 43, 82, 21, 61, 0, 66},
    //     {67, 42, 97, 52, 65, 88, 47, 88, 48, 76, 17, 66, 0}
    // };


    // // 14x14 Matrix
    // int matrix[N][N] = {
    //     {0, 67, 32, 89, 14, 76, 43, 91, 28, 55, 82, 37, 64, 95},
    //     {67, 0, 51, 23, 88, 41, 79, 16, 62, 94, 29, 73, 45, 18},
    //     {32, 51, 0, 74, 36, 92, 27, 85, 49, 11, 66, 53, 87, 22},
    //     {89, 23, 74, 0, 68, 15, 97, 44, 81, 33, 58, 90, 19, 71},
    //     {14, 88, 36, 68, 0, 59, 24, 77, 42, 86, 31, 63, 98, 26},
    //     {76, 41, 92, 15, 59, 0, 83, 38, 70, 25, 91, 47, 12, 84},
    //     {43, 79, 27, 97, 24, 83, 0, 61, 35, 88, 17, 72, 54, 29},
    //     {91, 16, 85, 44, 77, 38, 61, 0, 93, 52, 69, 21, 75, 48},
    //     {28, 62, 49, 81, 42, 70, 35, 93, 0, 64, 96, 39, 13, 87},
    //     {55, 94, 11, 33, 86, 25, 88, 52, 64, 0, 78, 46, 92, 31},
    //     {82, 29, 66, 58, 31, 91, 17, 69, 96, 78, 0, 84, 23, 65},
    //     {37, 73, 53, 90, 63, 47, 72, 21, 39, 46, 84, 0, 68, 15},
    //     {64, 45, 87, 19, 98, 12, 54, 75, 13, 92, 23, 68, 0, 77},
    //     {95, 18, 22, 71, 26, 84, 29, 48, 87, 31, 65, 15, 77, 0}
    // };


    // // 15x15 Matrix
    // int matrix[N][N] = {
    //     {0, 48, 76, 23, 91, 35, 82, 17, 64, 29, 87, 42, 95, 31, 53},
    //     {48, 0, 19, 85, 37, 92, 26, 71, 44, 98, 15, 67, 33, 79, 22},
    //     {76, 19, 0, 62, 88, 41, 96, 28, 73, 34, 81, 55, 12, 66, 47},
    //     {23, 85, 62, 0, 49, 77, 32, 94, 16, 83, 38, 90, 25, 58, 71},
    //     {91, 37, 88, 49, 0, 63, 21, 75, 52, 97, 29, 84, 46, 13, 68},
    //     {35, 92, 41, 77, 63, 0, 86, 24, 79, 11, 65, 39, 93, 27, 82},
    //     {82, 26, 96, 32, 21, 86, 0, 57, 45, 88, 14, 70, 36, 91, 19},
    //     {17, 71, 28, 94, 75, 24, 57, 0, 89, 33, 76, 22, 64, 48, 95},
    //     {64, 44, 73, 16, 52, 79, 45, 89, 0, 61, 27, 92, 18, 85, 31},
    //     {29, 98, 34, 83, 97, 11, 88, 33, 61, 0, 74, 56, 99, 42, 66},
    //     {87, 15, 81, 38, 29, 65, 14, 76, 27, 74, 0, 69, 23, 88, 51},
    //     {42, 67, 55, 90, 84, 39, 70, 22, 92, 56, 69, 0, 77, 35, 94},
    //     {95, 33, 12, 25, 46, 93, 36, 64, 18, 99, 23, 77, 0, 62, 28},
    //     {31, 79, 66, 58, 13, 27, 91, 48, 85, 42, 88, 35, 62, 0, 73},
    //     {53, 22, 47, 71, 68, 82, 19, 95, 31, 66, 51, 94, 28, 73, 0}
    // };




    TSPBnB<N> tsp(matrix);
    tsp.solve();
    tsp.printSolution();
    return 0;
}