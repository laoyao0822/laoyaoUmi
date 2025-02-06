cmake_minimum_required(VERSION 3.25)
project(laoyaoUmi)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_COMPILER "icx")
set(CMAKE_CXX_COMPILER "icpx")

# # 查找 htslib
# find_package(PkgConfig REQUIRED)
# pkg_check_modules(HTSLIB REQUIRED htslib)

# 指定库文件的路径
link_directories("/root/ASC25/m5C/htslib-intel_opt")

# 添加包含头文件的路径
include_directories("/root/ASC25/m5C/htslib-intel_opt/include")

# 查找 OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
    message("OpenMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found.")
endif()
message("${CMAKE_CXX_COMPILER_ID}")
# 针对 GCC 编译器配置 OpenMP 和优化
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    # 只添加一次 OpenMP 标志
    if (NOT CMAKE_CXX_FLAGS MATCHES "-fopenmp")
        add_compile_options(-qopenmp)
        message("GCC/Clang detected: Using --fopenmp")
    endif()
#-Ofast -mavx -mavx2 -msse -msse2
    # 设置 GCC 编译器优化等级为 -O3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fast")
    message("GCC/Clang detected: Using -O3 optimization")
else()
    message(WARNING "Unknown compiler, OpenMP flags might not be set.")
endif()


add_executable(
       laoyaoUmi main.cpp
)

# 链接到指定的 libhts.so
target_link_libraries(laoyaoUmi "/root/ASC25/m5C/htslib-intel_opt/libhts.so")

# 链接 htslib
target_include_directories(laoyaoUmi PRIVATE ${HTSLIB_INCLUDE_DIRS})
target_link_directories(laoyaoUmi PRIVATE ${HTSLIB_LIBRARY_DIRS})
target_link_libraries(laoyaoUmi PRIVATE ${HTSLIB_LIBRARIES})