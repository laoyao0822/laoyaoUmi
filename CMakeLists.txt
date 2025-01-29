cmake_minimum_required(VERSION 3.25)
project(laoyaoUmi)

set(CMAKE_CXX_STANDARD 17)

# 查找 htslib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

# 查找 OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_FOUND)
    message("OpenMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found.")
endif()

# 针对 GCC 编译器配置 OpenMP 和优化
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # 只添加一次 OpenMP 标志
    if (NOT CMAKE_CXX_FLAGS MATCHES "-fopenmp")
        add_compile_options(-fopenmp)
        message("GCC/Clang detected: Using -fopenmp")
    endif()

    # 设置 GCC 编译器优化等级为 -O3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    message("GCC/Clang detected: Using -O3 optimization")
else()
    message(WARNING "Unknown compiler, OpenMP flags might not be set.")
endif()


add_executable(
       laoyaoUmi main.cpp
)

# 链接 htslib
target_include_directories(laoyaoUmi PRIVATE ${HTSLIB_INCLUDE_DIRS})
target_link_directories(laoyaoUmi PRIVATE ${HTSLIB_LIBRARY_DIRS})
target_link_libraries(laoyaoUmi PRIVATE ${HTSLIB_LIBRARIES})