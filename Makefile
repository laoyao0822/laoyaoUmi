# 定义目标文件
TARGET = hisat-3n-table

# 定义源文件
SOURCE = hisat_3n_table.cpp

# 定义本地头文件
THREE_N_HEADERS = position_3n_table.h alignment_3n_table.h utility_3n_table.h

# 定义 htslib 的头文件路径
HTSLIB_INCLUDE = ../htslib-1.21/htslib

# 定义 htslib 的库文件路径
HTSLIB_LIB = ../htslib-1.21

# 定义编译器
CXX = /usr/bin/g++

# 定义编译器选项
CXXFLAGS = -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY -std=c++11 -DNDEBUG -fno-strict-aliasing -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -DBOWTIE_MM -pthread -I$(HTSLIB_INCLUDE)

# 定义链接器选项
LDFLAGS = -L$(HTSLIB_LIB) -lhts

# 定义其他宏定义
DEFINES = -DCOMPILER_OPTIONS="\"-O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY -std=c++11\"" \
          -DHISAT2_VERSION="\"`cat HISAT2_VERSION`\"" \
          -DBUILD_HOST="\"`hostname`\"" \
          -DBUILD_TIME="\"`date`\"" \
          -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\""

# 默认目标
all: $(TARGET)

# 编译目标
$(TARGET): $(SOURCE) $(THREE_N_HEADERS)
	$(CXX) $(CXXFLAGS) $(DEFINES) -o $@ $< $(LDFLAGS)

# 清理生成的文件
clean:
	rm -f $(TARGET)