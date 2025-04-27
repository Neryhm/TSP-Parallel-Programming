#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <cstring>
#include <omp.h>

class Node {
public:
    int level;
    std::vector<int> path;
    int** reducedMatrix;
    int cost;
    int lowerBound;
    int size;

    Node(int l, std::vector<int> p, int** m, int c, int lb, int s)
        : level(l), path(p), cost(c), lowerBound(lb), size(s) {
        reducedMatrix = new int*[s];
        for (int i = 0; i < s; ++i) {
            reducedMatrix[i] = new int[s];
            memcpy(reducedMatrix[i], m[i], s * sizeof(int));
        }
    }

    ~Node() {
        for (int i = 0; i < size; ++i) {
            delete[] reducedMatrix[i];
        }
        delete[] reducedMatrix;
    }
};

class TSPBranchAndBound {
private:
    int n;
    int** originalMatrix;
    int bestCost;
    std::vector<int> bestPath;
    std::ofstream logFile; // Log file stream

    struct CompareNode {
        bool operator()(const Node* a, const Node* b) {
            return a->lowerBound > b->lowerBound;
        }
    };

    std::priority_queue<Node*, std::vector<Node*>, CompareNode> pq;

    // Calculate MST cost for unvisited cities using Prim's algorithm
    int calculateMST(int** matrix, int size, const std::vector<int>& visited) {
        std::vector<bool> inMST(size, false);
        std::vector<int> key(size, INT_MAX);
        int mstCost = 0;

        for (int v : visited) {
            inMST[v] = true;
            key[v] = INT_MAX;
        }

        int start = 0;
        while (start < size && inMST[start]) {
            ++start;
        }
        if (start >= size) return 0;

        key[start] = 0;

        for (int count = 0; count < size - visited.size() - 1; ++count) {
            int minKey = INT_MAX, u = -1;
            for (int v = 0; v < size; ++v) {
                if (!inMST[v] && key[v] < minKey) {
                    minKey = key[v];
                    u = v;
                }
            }

            if (u == -1) break;

            inMST[u] = true;
            mstCost += minKey;

            for (int v = 0; v < size; ++v) {
                if (!inMST[v] && matrix[u][v] != INT_MAX && matrix[u][v] < key[v]) {
                    key[v] = matrix[u][v];
                }
            }
        }

        return mstCost;
    }

    int reduceMatrix(int** matrix, int size, const std::vector<int>& path) {
        int reduction = 0;

        for (int i = 0; i < size; ++i) {
            int minVal = INT_MAX;
            for (int j = 0; j < size; ++j) {
                if (matrix[i][j] < minVal) {
                    minVal = matrix[i][j];
                }
            }
            if (minVal != INT_MAX && minVal != 0) {
                reduction += minVal;
                for (int j = 0; j < size; ++j) {
                    if (matrix[i][j] != INT_MAX) {
                        matrix[i][j] -= minVal;
                    }
                }
            }
        }

        for (int j = 0; j < size; ++j) {
            int minVal = INT_MAX;
            for (int i = 0; i < size; ++i) {
                if (matrix[i][j] < minVal) {
                    minVal = matrix[i][j];
                }
            }
            if (minVal != INT_MAX && minVal != 0) {
                reduction += minVal;
                for (int i = 0; i < size; ++i) {
                    if (matrix[i][j] != INT_MAX) {
                        matrix[i][j] -= minVal;
                    }
                }
            }
        }

        if (path.size() < size) {
            reduction += calculateMST(matrix, size, path);
        }

        return reduction;
    }

    int** createNewMatrix(int parentCity, int childCity, int** parentMatrix) {
        int** newMatrix = new int*[n];
        for (int i = 0; i < n; ++i) {
            newMatrix[i] = new int[n];
            memcpy(newMatrix[i], parentMatrix[i], n * sizeof(int));
        }

        for (int j = 0; j < n; ++j) {
            newMatrix[parentCity][j] = INT_MAX;
        }

        for (int i = 0; i < n; ++i) {
            newMatrix[i][childCity] = INT_MAX;
        }

        newMatrix[childCity][0] = INT_MAX;

        return newMatrix;
    }

    void initializeRootNode() {
        int** rootMatrix = new int*[n];
        for (int i = 0; i < n; ++i) {
            rootMatrix[i] = new int[n];
            memcpy(rootMatrix[i], originalMatrix[i], n * sizeof(int));
        }

        std::vector<int> path = {0};
        int reduction = reduceMatrix(rootMatrix, n, path);
        Node* root = new Node(0, path, rootMatrix, 0, reduction, n);
        pq.push(root);

        for (int i = 0; i < n; ++i) {
            delete[] rootMatrix[i];
        }
        delete[] rootMatrix;

        std::cout << "Root node created. Lower Bound: " << reduction << std::endl;
        logFile << "Root node created. Lower Bound: " << reduction << std::endl;
    }

public:
    TSPBranchAndBound(int** matrix, int size) : n(size), bestCost(INT_MAX) {
        logFile.open("tsp_log_parallel.txt");
        if (!logFile.is_open()) {
            std::cerr << "Error opening log file!" << std::endl;
        }
        originalMatrix = new int*[n];
        for (int i = 0; i < n; ++i) {
            originalMatrix[i] = new int[n];
            memcpy(originalMatrix[i], matrix[i], n * sizeof(int));
        }
        initializeRootNode();
    }

    ~TSPBranchAndBound() {
        for (int i = 0; i < n; ++i) {
            delete[] originalMatrix[i];
        }
        delete[] originalMatrix;
        logFile.close();
    }

    void solve() {
        omp_set_dynamic(0);
        omp_set_num_threads(omp_get_max_threads());

        while (!pq.empty()) {
            Node* current = pq.top();
            pq.pop();

            if (current->lowerBound >= bestCost) {
                delete current;
                continue;
            }

            if (current->level == n - 1) {
                int finalCost = current->cost + originalMatrix[current->path.back()][0];
                if (finalCost < bestCost) {
                    #pragma omp critical
                    {
                        bestCost = finalCost;
                        bestPath = current->path;
                        bestPath.push_back(0);
                        std::cout << "New best path: ";
                        logFile << "New best path: ";
                        for (int city : bestPath) {
                            std::cout << city << " ";
                            logFile << city << " ";
                        }
                        std::cout << "| Cost: " << bestCost << std::endl;
                        logFile << "| Cost: " << bestCost << std::endl;
                    }
                }
                delete current;
                continue;
            }

            std::vector<int> candidates;
            for (int nextCity = 0; nextCity < n; ++nextCity) {
                if (std::find(current->path.begin(), current->path.end(), nextCity) == current->path.end()) {
                    candidates.push_back(nextCity);
                }
            }

            #pragma omp parallel for schedule(dynamic)
            for (int idx = 0; idx < (int)candidates.size(); ++idx) {
                int nextCity = candidates[idx];
                int** childMatrix = createNewMatrix(current->path.back(), nextCity, current->reducedMatrix);
                int edgeCost = originalMatrix[current->path.back()][nextCity];
                int newCost = current->cost + edgeCost;

                std::vector<int> newPath = current->path;
                newPath.push_back(nextCity);

                int reduction = reduceMatrix(childMatrix, n, newPath);
                int newBound = newCost + reduction;

                if (newBound < bestCost) {
                    #pragma omp critical
                    {
                        Node* child = new Node(current->level + 1, newPath, childMatrix, newCost, newBound, n);
                        pq.push(child);

                        std::cout << "Path: ";
                        logFile << "Path: ";
                        for (int city : newPath) {
                            std::cout << city << " ";
                            logFile << city << " ";
                        }
                        std::cout << "| Lower Bound: " << newBound << std::endl;
                        logFile << "| Lower Bound: " << newBound << std::endl;
                    }
                } else {
                    for (int i = 0; i < n; ++i) {
                        delete[] childMatrix[i];
                    }
                    delete[] childMatrix;
                }
            }
            delete current;
        }
    }

    int getBestCost() const { return bestCost; }
    const std::vector<int>& getBestPath() const { return bestPath; }

    void logFinalOutput(double duration) {
        std::cout << "\nOptimal Path: ";
        logFile << "\nOptimal Path: ";
        for (int city : getBestPath()) {
            std::cout << city << " ";
            logFile << city << " ";
        }
        std::cout << "\nTotal Cost: " << getBestCost() << std::endl;
        logFile << "\nTotal Cost: " << getBestCost() << std::endl;
        std::cout << "Runtime: " << duration << " seconds" << std::endl;
        logFile << "Runtime: " << duration << " seconds" << std::endl;
    }
};

int main() {
    const int n = 12;
    int matrix[n][n] = {
        {INT_MAX, 15, 18, 23, 13, 21, 17, 12, 20, 16, 24, 11},
        {14, INT_MAX, 19, 16, 22, 25, 11, 18, 24, 15, 21, 13},
        {17, 13, INT_MAX, 24, 15, 18, 23, 16, 21, 14, 22, 12},
        {21, 12, 19, INT_MAX, 17, 14, 26, 18, 22, 13, 20, 15},
        {23, 16, 20, 15, INT_MAX, 22, 13, 19, 17, 21, 12, 24},
        {18, 24, 17, 21, 12, INT_MAX, 15, 23, 14, 20, 16, 11},
        {13, 19, 22, 16, 25, 14, INT_MAX, 17, 21, 11, 23, 18},
        {20, 15, 18, 23, 17, 21, 12, INT_MAX, 16, 22, 13, 19},
        {22, 17, 21, 19, 14, 16, 24, 18, INT_MAX, 15, 20, 12},
        {16, 21, 14, 20, 22, 13, 17, 25, 12, INT_MAX, 19, 23},
        {24, 20, 23, 18, 15, 19, 21, 14, 17, 22, INT_MAX, 16},
        {11, 16, 12, 22, 19, 15, 20, 24, 13, 18, 21, INT_MAX}
    };

    int** dynamicMatrix = new int*[n];
    for (int i = 0; i < n; ++i) {
        dynamicMatrix[i] = new int[n];
        for (int j = 0; j < n; ++j) {
            dynamicMatrix[i][j] = matrix[i][j];
        }
    }

    TSPBranchAndBound tsp(dynamicMatrix, n);
    std::clock_t start = std::clock();
    tsp.solve();
    double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

    tsp.logFinalOutput(duration);

    for (int i = 0; i < n; ++i) {
        delete[] dynamicMatrix[i];
    }
    delete[] dynamicMatrix;

    return 0;
}