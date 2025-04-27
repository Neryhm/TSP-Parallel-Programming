#include <iostream>
#include <vector>
#include <queue>
#include <numeric>
#include <limits>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <chrono> // Include for timing

const int INF = std::numeric_limits<int>::max();
int N; // Number of cities
std::ofstream log_file;

// Structure to represent a node in the search tree
struct Node {
    std::vector<std::vector<int>> reduced_matrix;
    std::vector<int> path;
    int cost; // Lower bound
    int level; // Number of cities visited (level in the search tree)
    int city_id; // Current city ID

    // Overload the less than operator for priority queue (min-heap)
    bool operator>(const Node& other) const {
        return cost > other.cost;
    }
};

// Function to print the matrix (for debugging/logging)
void printMatrix(const std::vector<std::vector<int>>& matrix, const std::string& title = "Matrix") {
    log_file << title << ":\n";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (matrix[i][j] == INF) {
                log_file << std::setw(5) << "INF";
            } else {
                log_file << std::setw(5) << matrix[i][j];
            }
        }
        log_file << "\n";
    }
    log_file << std::endl;
}

// Function to reduce the cost matrix and return the reduction cost
int reduceMatrix(std::vector<std::vector<int>>& matrix) {
    int reduction_cost = 0;

    // Row reduction
    for (int i = 0; i < N; ++i) {
        int min_val = INF;
        for (int j = 0; j < N; ++j) {
            if (matrix[i][j] < min_val) {
                min_val = matrix[i][j];
            }
        }
        if (min_val != 0 && min_val != INF) {
            reduction_cost += min_val;
            for (int j = 0; j < N; ++j) {
                if (matrix[i][j] != INF) {
                    matrix[i][j] -= min_val;
                }
            }
        }
    }

    // Column reduction
    for (int j = 0; j < N; ++j) {
        int min_val = INF;
        for (int i = 0; i < N; ++i) {
            if (matrix[i][j] < min_val) {
                min_val = matrix[i][j];
            }
        }
        if (min_val != 0 && min_val != INF) {
            reduction_cost += min_val;
            for (int i = 0; i < N; ++i) {
                if (matrix[i][j] != INF) {
                    matrix[i][j] -= min_val;
                }
            }
        }
    }
    return reduction_cost;
}

// Function to create a new node
Node* createNode(std::vector<std::vector<int>> parent_matrix, const std::vector<int>& path,
                 int level, int current_city, int next_city) {
    Node* node = new Node;
    node->path = path;
    node->path.push_back(next_city);
    node->level = level;
    node->city_id = next_city;
    node->reduced_matrix = parent_matrix;

    // Set row of current_city and column of next_city to INF
    if (level > 0) { // Only if not the root node's initial setup
        for (int k = 0; k < N; ++k) {
            node->reduced_matrix[current_city][k] = INF; // Row
            node->reduced_matrix[k][next_city] = INF;   // Column
        }
        // Set matrix[next_city][start_city] to INF to prevent immediate return
        node->reduced_matrix[next_city][path[0]] = INF;
    }

    return node;
}

// Function to solve the TSP using Branch and Bound
void solveTSP(const std::vector<std::vector<int>>& initial_matrix) {
    N = initial_matrix.size();
    int min_cost = INF;
    std::vector<int> final_path;
    long long nodes_explored_count = 0; // Counter for loading indicator
    const int loading_update_frequency = 1000; // Print dot every 1000 nodes

    // Priority queue to store live nodes (min-heap based on cost)
    std::priority_queue<Node*, std::vector<Node*>, std::greater<Node*>> pq;

    // Create the root node
    Node* root = new Node;
    root->path.push_back(0); // Start at city 0
    root->level = 0;
    root->city_id = 0;
    root->reduced_matrix = initial_matrix;
    root->cost = reduceMatrix(root->reduced_matrix); // Initial lower bound

    log_file << "Initial Matrix Reduction:\n";
    log_file << "Initial Lower Bound (Cost): " << root->cost << "\n\n";

    pq.push(root);

    std::cout << "Solving TSP"; // Initial message for loading bar

    while (!pq.empty()) {
        // Get the node with the lowest cost
        Node* current_node = pq.top();
        pq.pop();
        nodes_explored_count++;

        // Loading indicator update
        if (nodes_explored_count % loading_update_frequency == 0) {
            std::cout << "." << std::flush;
        }

        // Simplified log for exploring a node
        log_file << "Explore: Level " << current_node->level << ", City " << current_node->city_id << ", Cost " << current_node->cost << ", Path: ";
        for (size_t i = 0; i < current_node->path.size(); ++i) log_file << current_node->path[i] << (i == current_node->path.size() - 1 ? "" : "-");
        log_file << "\n";

        // If the current node's cost is already greater than the minimum cost found so far, prune it
        if (current_node->cost >= min_cost) {
            log_file << "  Pruned (Cost " << current_node->cost << " >= Min " << min_cost << ")\n";
            delete current_node; // Free memory
            continue;
        }

        // If all cities have been visited (level == N - 1)
        if (current_node->level == N - 1) {
            // Complete the tour by returning to the start city
            current_node->path.push_back(current_node->path[0]); // Add start city to end
            log_file << "  Found Tour (Cost " << current_node->cost << "): ";
            for (size_t i = 0; i < current_node->path.size(); ++i) log_file << current_node->path[i] << (i == current_node->path.size() - 1 ? "" : " -> ");
            log_file << "\n";

            // Update the minimum cost if this tour is better
            if (current_node->cost < min_cost) {
                int old_min_cost = min_cost;
                min_cost = current_node->cost;
                final_path = current_node->path;
                log_file << "    *** New Best Tour! Cost updated from " << (old_min_cost == INF ? "INF" : std::to_string(old_min_cost)) << " to " << min_cost << " ***\n";
            }
            delete current_node; // Free memory
            continue; // Move to the next node in pq
        }

        // Branch: Explore child nodes for unvisited cities
        for (int next_city = 0; next_city < N; ++next_city) {
            // Check if the edge exists and the city hasn't been visited
            bool visited = false;
            for (int p : current_node->path) {
                if (p == next_city) {
                    visited = true;
                    break;
                }
            }

            if (current_node->reduced_matrix[current_node->city_id][next_city] != INF && !visited) {
                // Create child node
                Node* child_node = createNode(current_node->reduced_matrix, current_node->path,
                                              current_node->level + 1, current_node->city_id, next_city);

                // Calculate the cost for the child node
                int edge_cost = current_node->reduced_matrix[current_node->city_id][next_city];
                int reduction_cost = reduceMatrix(child_node->reduced_matrix);
                child_node->cost = current_node->cost + edge_cost + reduction_cost;

                // If the child's cost is less than the current minimum, add it to the queue
                if (child_node->cost < min_cost) {
                    pq.push(child_node);
                    log_file << "  Queueing Child: " << current_node->city_id << "->" << next_city << " (Cost " << child_node->cost << ")\n";
                } else {
                    log_file << "  Pruning Child: " << current_node->city_id << "->" << next_city << " (Cost " << child_node->cost << " >= Min " << min_cost << ")\n";
                    delete child_node; // Free memory
                }
            }
        }
        delete current_node; // Free memory of the processed node
    }
    std::cout << std::endl; // Newline after loading indicator finishes

    // Print the final result
    std::cout << "Nodes Explored: " << nodes_explored_count << std::endl;
    std::cout << "Minimum Cost: " << min_cost << std::endl;
    std::cout << "Optimal Path: ";
    for (int i = 0; i < final_path.size(); ++i) {
        std::cout << final_path[i] << (i == final_path.size() - 1 ? "" : " -> ");
    }
    std::cout << std::endl;

    log_file << "========================================\n";
    log_file << "           Final Result\n";
    log_file << "========================================\n";
    log_file << "Nodes Explored: " << nodes_explored_count << std::endl;
    log_file << "Minimum Cost: " << min_cost << std::endl;
    log_file << "Optimal Path: ";
    for (int i = 0; i < final_path.size(); ++i) {
        log_file << final_path[i] << (i == final_path.size() - 1 ? "" : " -> ");
    }
    log_file << std::endl;
}

int main() {
    // Example 19x19 symmetric cost matrix
    std::vector<std::vector<int>> cost_matrix = {
        {INF, 15, 18, 23, 13, 21, 17, 12, 20, 16, 24, 11, 19, 14, 22, 25, 10, 16, 18},
        {15, INF, 19, 16, 22, 25, 11, 18, 24, 15, 21, 13, 17, 20, 12, 14, 23, 19, 15},
        {18, 19, INF, 24, 15, 18, 23, 16, 21, 14, 22, 12, 20, 13, 17, 19, 11, 24, 16},
        {23, 16, 24, INF, 17, 14, 26, 18, 22, 13, 20, 15, 12, 21, 19, 23, 14, 17, 25},
        {13, 22, 15, 17, INF, 22, 13, 19, 17, 21, 12, 24, 16, 18, 20, 11, 25, 14, 23},
        {21, 25, 18, 14, 22, INF, 15, 23, 14, 20, 16, 11, 19, 12, 24, 17, 13, 21, 18},
        {17, 11, 23, 26, 13, 15, INF, 17, 21, 11, 23, 18, 14, 20, 16, 22, 19, 12, 24},
        {12, 18, 16, 18, 19, 23, 17, INF, 16, 22, 13, 19, 21, 15, 11, 24, 20, 14, 17},
        {20, 24, 21, 22, 17, 14, 21, 16, INF, 15, 20, 12, 18, 23, 19, 13, 11, 25, 22},
        {16, 15, 14, 13, 21, 20, 11, 22, 15, INF, 19, 23, 17, 12, 24, 18, 20, 16, 21},
        {24, 21, 22, 20, 12, 16, 23, 13, 20, 19, INF, 16, 14, 18, 15, 21, 17, 11, 23},
        {11, 13, 12, 15, 24, 11, 18, 19, 12, 23, 16, INF, 20, 22, 17, 14, 25, 19, 13},
        {19, 17, 20, 12, 16, 19, 14, 21, 18, 17, 14, 20, INF, 15, 23, 11, 22, 24, 16},
        {14, 20, 13, 21, 18, 12, 20, 15, 23, 12, 18, 22, 15, INF, 19, 16, 13, 21, 17},
        {22, 12, 17, 19, 20, 24, 16, 11, 19, 24, 15, 17, 23, 19, INF, 18, 14, 20, 12},
        {25, 14, 19, 23, 11, 17, 22, 24, 13, 18, 21, 14, 11, 16, 18, INF, 20, 15, 19},
        {10, 23, 11, 14, 25, 13, 19, 20, 11, 20, 17, 25, 22, 13, 14, 20, INF, 18, 16},
        {16, 19, 24, 17, 14, 21, 12, 14, 25, 16, 11, 19, 24, 21, 20, 15, 18, INF, 22},
        {18, 15, 16, 25, 23, 18, 24, 17, 22, 21, 23, 13, 16, 17, 12, 19, 16, 22, INF}
    };

    // Open log file
    log_file.open("tsp_log_sequential.txt");
    if (!log_file.is_open()) {
        std::cerr << "Error opening log file!" << std::endl;
        return 1;
    }

    log_file << "Starting TSP Solver (Branch and Bound with Matrix Reduction)\n";
    log_file << "Number of cities: " << cost_matrix.size() << "\n";
    printMatrix(cost_matrix, "Initial Cost Matrix");

    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();

    solveTSP(cost_matrix);

    // Stop timer
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // Log runtime
    log_file << "========================================\n";
    log_file << "           Performance\n";
    log_file << "========================================\n";
    log_file << "Total Runtime: " << duration.count() << " ms" << std::endl;
    std::cout << "Total Runtime: " << duration.count() << " ms" << std::endl;

    log_file.close();

    return 0;
}
