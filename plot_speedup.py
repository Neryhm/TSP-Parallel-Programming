import pandas as pd
import matplotlib
# Đặt chế độ sử dụng matplotlib để không cần hiển thị cửa sổ đồ họa
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np

# Đọc dữ liệu từ file CSV (nếu bạn tạo file CSV)
# Giả sử file CSV có các cột: 'Processors_OpenMP', 'Speedup_OpenMP'
# và một điểm dữ liệu riêng cho CUDA.

# --- Cách 1: Nhập dữ liệu thủ công (nếu ít điểm dữ liệu) ---
processors_openmp = np.array([1, 2, 4, 8]) # Số lượng bộ xử lý/luồng cho OpenMP
# Thời gian chạy tuần tự (hoặc OpenMP với 1 luồng) - ví dụ: 34600 ms
time_sequential = 43250.0

# Thời gian chạy OpenMP tương ứng với processors_openmp - ví dụ
times_openmp = np.array([34946.0, 21064.0, 19314.0, 17232.0])

# Thời gian chạy CUDA - ví dụ (coi như một điểm)
time_cuda = 7211.0 # Thời gian trung bình bạn đo được
processors_cuda_point = 8 # Vị trí để đánh dấu điểm CUDA trên trục x (có thể điều chỉnh)

# Tính Speedup cho OpenMP
speedup_openmp = time_sequential / times_openmp

# Tính Speedup cho CUDA (so với tuần tự)
speedup_cuda = time_sequential / time_cuda

# --- Cách 2: Đọc từ file CSV (khuyến khích nếu nhiều dữ liệu) ---
# # Tạo file 'results.csv' với nội dung ví dụ:
# # Processors,Time_OpenMP
# # 1,34600
# # 2,18000
# # 4,9000
# # 8,4800
# # 16,3460
# # 32,2500
# # 64,2000
# #
# # # Dữ liệu CUDA (có thể thêm vào file khác hoặc xử lý riêng)
# # # CUDA_Time,CUDA_Processors_Equivalent (ví dụ)
# # # 2000,64

# try:
#     data_openmp = pd.read_csv('results_openmp.csv') # Giả sử có file riêng cho OpenMP
#     processors_openmp = data_openmp['Processors'].values
#     times_openmp = data_openmp['Time_OpenMP'].values
#     time_sequential_from_file = times_openmp[0] # Giả sử dòng đầu tiên là tuần tự (1 processor)
#     speedup_openmp = time_sequential_from_file / times_openmp

#     # Đọc dữ liệu CUDA nếu có file riêng hoặc thêm cột vào file chính
#     # time_cuda = 2000.0 # Lấy từ file hoặc nhập tay
#     # processors_cuda_point = 64
#     # speedup_cuda = time_sequential_from_file / time_cuda

# except FileNotFoundError:
#     print("Không tìm thấy file CSV, sử dụng dữ liệu mẫu thủ công.")
#     # Dùng lại dữ liệu thủ công từ Cách 1 nếu file không tồn tại


# --- Vẽ biểu đồ ---
plt.figure(figsize=(10, 6))

# Biểu đồ Speedup thực tế (OpenMP)
plt.plot(processors_openmp, speedup_openmp, marker='o', linestyle='-', color='magenta', label='Actual Speedup (OpenMP)')

# Đánh dấu điểm CUDA
# Bạn có thể chọn cách thể hiện điểm CUDA:
# 1. Là một điểm riêng biệt nếu bạn coi nó là một cấu hình cụ thể.
# 2. Nếu bạn có nhiều điểm dữ liệu CUDA (ví dụ: với các BATCH_SIZE khác nhau), bạn có thể vẽ một đường riêng.
# Ở đây, ta vẽ một điểm duy nhất.
plt.plot(processors_cuda_point, speedup_cuda, marker='s', color='blue', markersize=8, linestyle='None', label=f'CUDA Speedup (at {processors_cuda_point} "equiv." proc.)')
plt.text(processors_cuda_point, speedup_cuda + 0.5, f'CUDA: {speedup_cuda:.1f}x', color='blue')


# Biểu đồ Speedup lý tưởng (y=x)
# Chọn giá trị lớn nhất của processors để vẽ đường lý tưởng
max_processors = np.max(processors_openmp)
if processors_cuda_point > max_processors:
    max_processors = processors_cuda_point

ideal_processors = np.array([1, max_processors if max_processors > 1 else 2]) # Cần ít nhất 2 điểm để vẽ đường thẳng
ideal_speedup = ideal_processors
plt.plot(ideal_processors, ideal_speedup, linestyle='--', color='blue', label='Ideal Speedup (y=x)')


# Cài đặt biểu đồ
plt.title('Speedup vs Number of Processors/Threads')
plt.xlabel('Number of Processors / Threads')
plt.ylabel('Speedup')
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5) # Thêm lưới cho cả trục chính và phụ (nếu có log scale)

# Thiết lập giới hạn cho trục (tùy chọn, điều chỉnh cho phù hợp với dữ liệu của bạn)
plt.xlim(left=0.8, right=max_processors + 5 if max_processors > 1 else 2.2) # Để điểm 1 không bị sát lề
min_speedup_display = 0.8
max_speedup_display = max(np.max(speedup_openmp) if len(speedup_openmp) > 0 else 0, speedup_cuda if speedup_cuda else 0, max_processors if max_processors > 1 else 2)
plt.ylim(bottom=min_speedup_display, top=max_speedup_display + 2)


# --- Tùy chọn: Trục Logarit (giống biểu đồ thứ 2 trong ảnh của bạn) ---
# Nếu muốn sử dụng trục logarit cho cả x và y:
# plt.xscale('log')
# plt.yscale('log')
# # Khi dùng log scale, cần đảm bảo processors không bắt đầu từ 0.
# # Giới hạn trục và các điểm dữ liệu cần được điều chỉnh cho phù hợp.
# # Ví dụ, nếu dùng log scale:
# if np.min(processors_openmp) < 1: # Đảm bảo không có giá trị <=0 cho trục log
#      print("Warning: Log scale requires processor counts > 0.")
# else:
#      plt.xlim(left=np.min(processors_openmp) * 0.8, right=max_processors * 1.2)
#      plt.ylim(bottom=0.8, top=max_speedup_display * 1.2)


# Hiển thị biểu đồ
plt.show()

# Lưu biểu đồ ra file
plt.savefig('speedup_comparison_chart.png')
