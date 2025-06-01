#include <iostream>
#include <vector>
#include <queue>
#include <numeric>
#include <limits>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <chrono>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// ==================================================================================
// Constants and Global (Host) Variables
// ==================================================================================
const int INF = std::numeric_limits<int>::max();
const int MAX_CITIES = 20; // Giới hạn số thành phố, điều chỉnh nếu cần

int N_actual; // Số thành phố thực tế của bài toán
std::ofstream log_file_cuda;

// ==================================================================================
// CUDA Error Checking Utility
// ==================================================================================
#define cudaErrCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// ==================================================================================
// Data Structures
// ==================================================================================
struct Node {
    int reduced_matrix[MAX_CITIES * MAX_CITIES]; // Ma trận rút gọn (mảng 1D, stride MAX_CITIES)
    int path[MAX_CITIES];
    int cost;
    int level;
    int city_id;

    // Hàm khởi tạo cho cả host và device
    __host__ __device__ Node() : cost(0), level(0), city_id(0) {
        for(int i=0; i < MAX_CITIES; ++i) path[i] = -1;
        for(int i=0; i < MAX_CITIES * MAX_CITIES; ++i) reduced_matrix[i] = INF; // Khởi tạo ma trận bằng INF
    }

    // Để sử dụng với std::priority_queue trên host
    bool operator>(const Node& other) const {
        return cost > other.cost;
    }
};


// ==================================================================================
// Device Code (__device__ functions)
// ==================================================================================

// Hàm rút gọn ma trận trên device
__device__ int reduceMatrix_device(int* matrix, int n_val) { // n_val là N_actual
    int reduction_cost = 0;
    const int stride = MAX_CITIES; // LUÔN dùng MAX_CITIES làm stride

    // Row reduction
    for (int i = 0; i < n_val; ++i) { // Lặp qua N_actual hàng
        int min_val_row = INF;
        for (int j = 0; j < n_val; ++j) { // Lặp qua N_actual cột
            if (matrix[i * stride + j] < min_val_row) { // Dùng stride
                min_val_row = matrix[i * stride + j];
            }
        }
        if (min_val_row != 0 && min_val_row != INF) {
            reduction_cost += min_val_row;
            for (int j = 0; j < n_val; ++j) {
                if (matrix[i * stride + j] != INF) {
                    matrix[i * stride + j] -= min_val_row;
                }
            }
        }
    }

    // Column reduction
    for (int j = 0; j < n_val; ++j) { // Lặp qua N_actual cột
        int min_val_col = INF;
        for (int i = 0; i < n_val; ++i) { // Lặp qua N_actual hàng
            if (matrix[i * stride + j] < min_val_col) { // Dùng stride
                min_val_col = matrix[i * stride + j];
            }
        }
        if (min_val_col != 0 && min_val_col != INF) {
            reduction_cost += min_val_col;
            for (int i = 0; i < n_val; ++i) {
                if (matrix[i * stride + j] != INF) {
                    matrix[i * stride + j] -= min_val_col;
                }
            }
        }
    }
    return reduction_cost;
}

// ==================================================================================
// Kernel Code (__global__ functions)
// ==================================================================================
__global__ void processNodes_kernel(Node* d_input_nodes,
                                   int num_nodes_in_batch,
                                   int current_min_cost_host,
                                   int n_val,                   // N_actual
                                   Node* d_output_children,
                                   int* d_children_count,
                                   int max_children_per_node)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = MAX_CITIES; // Sử dụng MAX_CITIES làm stride

    if (idx < num_nodes_in_batch) {
        Node current_node = d_input_nodes[idx]; // Lấy một bản sao để làm việc

        if (current_node.cost >= current_min_cost_host) {
            return;
        }

        if (current_node.level == n_val - 1) {
             int output_idx = atomicAdd(d_children_count, 1);
             if (output_idx < num_nodes_in_batch * max_children_per_node) {
                d_output_children[output_idx] = current_node;
             }
            return;
        }

        for (int next_city = 0; next_city < n_val; ++next_city) {
            bool visited = false;
            for (int i = 0; i <= current_node.level; ++i) {
                if (current_node.path[i] == next_city) {
                    visited = true;
                    break;
                }
            }

            if (current_node.reduced_matrix[current_node.city_id * stride + next_city] != INF && !visited) {
                Node child_node; // Gọi hàm khởi tạo __device__ Node() -> ma trận được init bằng INF

                // Sao chép phần N_actual x N_actual của ma trận cha
                for (int r = 0; r < n_val; ++r) {
                    for (int c = 0; c < n_val; ++c) {
                        child_node.reduced_matrix[r * stride + c] = current_node.reduced_matrix[r * stride + c];
                    }
                }
                // Phần padding của child_node.reduced_matrix (nếu n_val < MAX_CITIES) đã được init là INF bởi constructor.

                for(int i=0; i <= current_node.level; ++i) child_node.path[i] = current_node.path[i];

                child_node.level = current_node.level + 1;
                child_node.path[child_node.level] = next_city;
                child_node.city_id = next_city;

                int edge_cost_in_parent_matrix = current_node.reduced_matrix[current_node.city_id * stride + next_city];

                // Đặt các hàng/cột thành INF trong ma trận con, dùng stride
                // Chỉ cần thao tác trên phần N_actual x N_actual của ma trận con
                for (int k = 0; k < n_val; ++k) {
                    child_node.reduced_matrix[current_node.city_id * stride + k] = INF; // Hàng của thành phố hiện tại
                    child_node.reduced_matrix[k * stride + next_city] = INF;           // Cột của thành phố kế tiếp
                }

                if (child_node.path[0] != -1) { // path[0] là thành phố bắt đầu của tour
                   child_node.reduced_matrix[next_city * stride + child_node.path[0]] = INF; // Ngăn quay lại ngay TP bắt đầu tour
                }


                int reduction_cost_child = reduceMatrix_device(child_node.reduced_matrix, n_val);
                child_node.cost = current_node.cost + edge_cost_in_parent_matrix + reduction_cost_child;

                if (child_node.cost < current_min_cost_host) {
                    int output_idx = atomicAdd(d_children_count, 1);
                    if (output_idx < num_nodes_in_batch * max_children_per_node) {
                        d_output_children[output_idx] = child_node;
                    }
                }
            }
        }
    }
}


// ==================================================================================
// Host Code
// ==================================================================================

void printMatrix_host(const int* matrix_flat, int n_val_to_print, const std::string& title = "Matrix") {
    log_file_cuda << title << " (Printing " << n_val_to_print << "x" << n_val_to_print << " view, Stored Stride=" << MAX_CITIES << "):\n";
    for (int i = 0; i < n_val_to_print; ++i) {
        for (int j = 0; j < n_val_to_print; ++j) {
            if (matrix_flat[i * MAX_CITIES + j] == INF) { // Dùng MAX_CITIES làm stride
                log_file_cuda << std::setw(5) << "INF";
            } else {
                log_file_cuda << std::setw(5) << matrix_flat[i * MAX_CITIES + j];
            }
        }
        log_file_cuda << "\n";
    }
    log_file_cuda << std::endl;
}

void printMatrix_host_vec(const std::vector<std::vector<int>>& matrix_vec, int n_val, const std::string& title = "Matrix") {
    log_file_cuda << title << " (N=" << n_val << "):\n";
    for (int i = 0; i < n_val; ++i) {
        for (int j = 0; j < n_val; ++j) {
            if (matrix_vec[i][j] == INF) {
                log_file_cuda << std::setw(5) << "INF";
            } else {
                log_file_cuda << std::setw(5) << matrix_vec[i][j];
            }
        }
        log_file_cuda << "\n";
    }
    log_file_cuda << std::endl;
}

int reduceMatrix_host_vec(std::vector<std::vector<int>>& matrix, int n_val) {
    int reduction_cost = 0;
    for (int i = 0; i < n_val; ++i) {
        int min_val = INF;
        for (int j = 0; j < n_val; ++j) if (matrix[i][j] < min_val) min_val = matrix[i][j];
        if (min_val != 0 && min_val != INF) {
            reduction_cost += min_val;
            for (int j = 0; j < n_val; ++j) if (matrix[i][j] != INF) matrix[i][j] -= min_val;
        }
    }
    for (int j = 0; j < n_val; ++j) {
        int min_val = INF;
        for (int i = 0; i < n_val; ++i) if (matrix[i][j] < min_val) min_val = matrix[i][j];
        if (min_val != 0 && min_val != INF) {
            reduction_cost += min_val;
            for (int i = 0; i < n_val; ++i) if (matrix[i][j] != INF) matrix[i][j] -= min_val;
        }
    }
    return reduction_cost;
}

void solveTSP_cuda(const std::vector<std::vector<int>>& initial_matrix_vec) {
    N_actual = initial_matrix_vec.size();
    if (N_actual > MAX_CITIES) {
        std::cerr << "Error: N_actual (" << N_actual << ") exceeds MAX_CITIES (" << MAX_CITIES << "). Please increase MAX_CITIES and recompile." << std::endl;
        return;
    }

    int min_cost_host = INF;
    std::vector<int> final_path_host(N_actual + 1);
    long long nodes_processed_on_gpu_count = 0; // Đếm số nút thực sự được gửi lên GPU và xử lý

    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq_host;

    Node root_host_node; // Constructor sẽ khởi tạo ma trận bằng INF
    root_host_node.level = 0;
    root_host_node.city_id = 0;
    root_host_node.path[0] = 0;

    std::vector<std::vector<int>> temp_matrix_for_reduction = initial_matrix_vec;
    root_host_node.cost = reduceMatrix_host_vec(temp_matrix_for_reduction, N_actual);

    // Sao chép phần N_actual x N_actual từ temp_matrix_for_reduction
    // vào root_host_node.reduced_matrix sử dụng MAX_CITIES làm stride
    for(int i=0; i < N_actual; ++i) {
        for(int j=0; j < N_actual; ++j) {
            root_host_node.reduced_matrix[i * MAX_CITIES + j] = temp_matrix_for_reduction[i][j];
        }
    }
    // Phần padding của root_host_node.reduced_matrix đã được init là INF bởi constructor.

    pq_host.push(root_host_node);
    log_file_cuda << "Initial Root Node Cost (Host): " << root_host_node.cost << std::endl;
    // printMatrix_host(root_host_node.reduced_matrix, N_actual, "Initial Root Matrix");


    const int BATCH_SIZE = 1024;
    std::vector<Node> batch_nodes_host(BATCH_SIZE);
    std::vector<Node> children_nodes_from_gpu_host;

    Node* d_input_nodes;
    Node* d_output_children;
    int* d_children_count;

    cudaErrCheck(cudaMalloc((void**)&d_input_nodes, BATCH_SIZE * sizeof(Node)));
    int max_children_estimate = BATCH_SIZE * (N_actual > 1 ? N_actual - 1 : 1);
    if (max_children_estimate == 0 && N_actual == 1) max_children_estimate = BATCH_SIZE;
    cudaErrCheck(cudaMalloc((void**)&d_output_children, max_children_estimate * sizeof(Node)));
    cudaErrCheck(cudaMalloc((void**)&d_children_count, sizeof(int)));

    std::cout << "Solving TSP with CUDA (N=" << N_actual << ", Batch Size=" << BATCH_SIZE << ")" << std::endl;

    while(!pq_host.empty()) {
        int current_batch_fill_size = 0;
        for(int i=0; i < BATCH_SIZE && !pq_host.empty(); ++i) {
            Node top_node = pq_host.top();
            pq_host.pop();

            if (top_node.cost >= min_cost_host) {
                i--; // Không tính nút này, thử lấy nút khác
                continue;
            }
            batch_nodes_host[i] = top_node;
            current_batch_fill_size++;
        }

        if (current_batch_fill_size == 0) break;

        nodes_processed_on_gpu_count += current_batch_fill_size;
        if (nodes_processed_on_gpu_count % 10000 < current_batch_fill_size) {
             std::cout << "." << std::flush;
        }

        cudaErrCheck(cudaMemcpy(d_input_nodes, batch_nodes_host.data(), current_batch_fill_size * sizeof(Node), cudaMemcpyHostToDevice));
        int h_children_count_reset = 0;
        cudaErrCheck(cudaMemcpy(d_children_count, &h_children_count_reset, sizeof(int), cudaMemcpyHostToDevice));

        int threads_per_block = 256;
        int blocks_per_grid = (current_batch_fill_size + threads_per_block - 1) / threads_per_block;

        processNodes_kernel<<<blocks_per_grid, threads_per_block>>>(
            d_input_nodes,
            current_batch_fill_size,
            min_cost_host,
            N_actual,
            d_output_children,
            d_children_count,
            (N_actual > 1 ? N_actual - 1 : 1) 
        );
        cudaErrCheck(cudaGetLastError());
        cudaErrCheck(cudaDeviceSynchronize());

        int h_children_count_from_gpu;
        cudaErrCheck(cudaMemcpy(&h_children_count_from_gpu, d_children_count, sizeof(int), cudaMemcpyDeviceToHost));

        if (h_children_count_from_gpu > 0) {
            // Đảm bảo children_nodes_from_gpu_host đủ lớn
            if (children_nodes_from_gpu_host.size() < (size_t)h_children_count_from_gpu) {
                 children_nodes_from_gpu_host.resize(h_children_count_from_gpu);
            }
            cudaErrCheck(cudaMemcpy(children_nodes_from_gpu_host.data(), d_output_children, h_children_count_from_gpu * sizeof(Node), cudaMemcpyDeviceToHost));

            for (int i=0; i < h_children_count_from_gpu; ++i) {
                Node child = children_nodes_from_gpu_host[i];
                if (child.cost >= min_cost_host) continue;

                if (child.level == N_actual - 1) {
                    int actual_tour_cost = child.cost;
                    if (initial_matrix_vec[child.city_id][child.path[0]] != INF) { // Kiểm tra cạnh cuối có tồn tại trong ma trận gốc
                        if (actual_tour_cost < min_cost_host) {
                            min_cost_host = actual_tour_cost;
                            for(int k=0; k <= child.level; ++k) final_path_host[k] = child.path[k];
                            final_path_host[N_actual] = child.path[0];
                            log_file_cuda << "  *** New Best Tour (GPU)! Cost: " << min_cost_host << " *** Path: ";
                            for(int k=0; k<=N_actual; ++k) log_file_cuda << final_path_host[k] << (k==N_actual ? "" : " -> ");
                            log_file_cuda << std::endl;
                        }
                    }
                } else {
                    pq_host.push(child);
                }
            }
        }
    }
    std::cout << std::endl << "Done." << std::endl;

    cudaErrCheck(cudaFree(d_input_nodes));
    cudaErrCheck(cudaFree(d_output_children));
    cudaErrCheck(cudaFree(d_children_count));

    std::cout << "Nodes processed on GPU (approx): " << nodes_processed_on_gpu_count << std::endl;
    std::cout << "Minimum Cost: " << min_cost_host << std::endl;
    std::cout << "Optimal Path: ";
    for (int i = 0; i <= N_actual; ++i) {
        std::cout << final_path_host[i] << (i == N_actual ? "" : " -> ");
    }
    std::cout << std::endl;

    log_file_cuda << "========================================\n";
    log_file_cuda << "           Final Result (CUDA)\n";
    log_file_cuda << "========================================\n";
    log_file_cuda << "Nodes processed on GPU (approx): " << nodes_processed_on_gpu_count << std::endl;
    log_file_cuda << "Minimum Cost: " << min_cost_host << std::endl;
    log_file_cuda << "Optimal Path: ";
    for (int i = 0; i <= N_actual; ++i) {
        log_file_cuda << final_path_host[i] << (i == N_actual ? "" : " -> ");
    }
    log_file_cuda << std::endl;
}

int main() {
    // std::vector<std::vector<int>> cost_matrix = {
    //     {INF, 10, 15, 20},
    //     {10, INF, 35, 25},
    //     {15, 35, INF, 30},
    //     {20, 25, 30, INF}
    // };
    // std::vector<std::vector<int>> cost_matrix = { // Test case 2 (Answer: 25)
    //     {INF, 20, 30, 10, 11},
    //     {15, INF, 16,  4,  2},
    //     { 3,  5, INF,  2,  4},
    //     {19,  6, 18, INF,  3},
    //     {16,  4,  7, 16, INF}
    // };
     std::vector<std::vector<int>> cost_matrix = { // Test case 3 from sequential (Answer: 242)
        {INF, 15, 18, 23, 13, 21, 17, 12, 20, 16, 24, 11, 19, 14, 22, 25, 10, 16, 18, 36},
        {15, INF, 19, 16, 22, 25, 11, 18, 24, 15, 21, 13, 17, 20, 12, 14, 23, 19, 15, 32},
        {18, 19, INF, 24, 15, 18, 23, 16, 21, 14, 22, 12, 20, 13, 17, 19, 11, 24, 16, 27},
        {23, 16, 24, INF, 17, 14, 26, 18, 22, 13, 20, 15, 12, 21, 19, 23, 14, 17, 25, 30},
        {13, 22, 15, 17, INF, 22, 13, 19, 17, 21, 12, 24, 16, 18, 20, 11, 25, 14, 23, 20},
        {21, 25, 18, 14, 22, INF, 15, 23, 14, 20, 16, 11, 19, 12, 24, 17, 13, 21, 18, 21},
        {17, 11, 23, 26, 13, 15, INF, 17, 21, 11, 23, 18, 14, 20, 16, 22, 19, 12, 24, 12},
        {12, 18, 16, 18, 19, 23, 17, INF, 16, 22, 13, 19, 21, 15, 11, 24, 20, 14, 17, 28},
        {20, 24, 21, 22, 17, 14, 21, 16, INF, 15, 20, 12, 18, 23, 19, 13, 11, 25, 22, 35},
        {16, 15, 14, 13, 21, 20, 11, 22, 15, INF, 19, 23, 17, 12, 24, 18, 20, 16, 21, 9},
        {24, 21, 22, 20, 12, 16, 23, 13, 20, 19, INF, 16, 14, 18, 15, 21, 17, 11, 23, 41},
        {11, 13, 12, 15, 24, 11, 18, 19, 12, 23, 16, INF, 20, 22, 17, 14, 25, 19, 13, 23},
        {19, 17, 20, 12, 16, 19, 14, 21, 18, 17, 14, 20, INF, 15, 23, 11, 22, 24, 16, 22},
        {14, 20, 13, 21, 18, 12, 20, 15, 23, 12, 18, 22, 15, INF, 19, 16, 13, 21, 17, 18},
        {22, 12, 17, 19, 20, 24, 16, 11, 19, 24, 15, 17, 23, 19, INF, 18, 14, 20, 12, 25},
        {25, 14, 19, 23, 11, 17, 22, 24, 13, 18, 21, 14, 11, 16, 18, INF, 20, 15, 19, 16},
        {10, 23, 11, 14, 25, 13, 19, 20, 11, 20, 17, 25, 22, 13, 14, 20, INF, 18, 16, 31},
        {16, 19, 24, 17, 14, 21, 12, 14, 25, 16, 11, 19, 24, 21, 20, 15, 18, INF, 22, 33},
        {18, 15, 16, 25, 23, 18, 24, 17, 22, 21, 23, 13, 16, 17, 12, 19, 16, 22, INF, 29},
        {36, 32, 27, 30, 20, 21, 12, 28, 35, 9, 41, 23, 22, 18, 25, 16, 31, 33, 29, INF}
     };


    if (cost_matrix.empty() || cost_matrix.size() > MAX_CITIES) {
        std::cerr << "Cost matrix is empty or too large for MAX_CITIES (Current MAX_CITIES = " << MAX_CITIES <<").\n";
        return 1;
    }

    log_file_cuda.open("tsp_log_cuda.txt");
    if (!log_file_cuda.is_open()) {
        std::cerr << "Error opening CUDA log file!" << std::endl;
        return 1;
    }

    log_file_cuda << "Starting TSP Solver (CUDA Hybrid)\n";
    log_file_cuda << "Number of cities: " << cost_matrix.size() << "\n";
    printMatrix_host_vec(cost_matrix, cost_matrix.size(), "Initial Cost Matrix (Host)");

    auto start_time = std::chrono::high_resolution_clock::now();
    solveTSP_cuda(cost_matrix);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Total Runtime (CUDA): " << duration.count() << " ms" << std::endl;
    log_file_cuda << "========================================\n";
    log_file_cuda << "           Performance\n";
    log_file_cuda << "========================================\n";
    log_file_cuda << "Total Runtime: " << duration.count() << " ms" << std::endl;

    log_file_cuda.close();

    return 0;
}