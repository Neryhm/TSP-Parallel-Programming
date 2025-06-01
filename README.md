For Sequential: g++ -fopenmp Sequential/TSP_Sequential.cpp -o Sequential/TSP_Sequential
For OpenMP: g++ -fopenmp OpenMP/TSP_Omp.cpp -o OpenMP/TSP_Omp

Parallelizing Step-by-step plan:
1. Include OpenMP Header: Add the <omp.h> header file.
2. Modify solveTSP Function:
- Wrap the main processing loop (while (!pq.empty())) within an OpenMP parallel region (#pragma omp parallel).
- Protect shared data access (pq, min_cost, final_path, log_file, nodes_explored_count) using OpenMP critical sections (#pragma omp critical) or atomic operations (#pragma omp atomic).
- Each thread will repeatedly attempt to pop a node from the shared priority queue, process it, and potentially add new child nodes back to the queue.
- The parallel region will implicitly handle thread management and termination when all threads are idle and the queue is empty.
- Remove the sequential loading indicator (.) as it's not suitable for parallel execution.
3. Modify main Function:
- Update the log file name to indicate OpenMP execution (e.g., tsp_log_openmp.txt).
- Optionally, print the maximum number of threads used by OpenMP.