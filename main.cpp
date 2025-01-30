//
// Created by laoyao on 2025/1/27.
//
#include "htslib/sam.h"
#include <iostream>
#include <omp.h>
#include <regex>
#include <bitset>
#include <string>
#include <stdexcept>
#include <unordered_map>
using namespace std;
int get_avg_qual(bam1_t* b);
// Initialize the static member
std::string extractUMI(const std::string& readName, char sep);
int get_unclipped_start(bam1_t *b) ;
int get_unclipped_end(bam1_t *b);

const char* get_reference_name(bam1_t *b, sam_hdr_t *header);
/**
 * 表示位置的结构体
 */
struct Alignment {
    bool isReversed;
    int pos;
    string refName;
    bool operator==(const Alignment &other) const {
        return isReversed == other.isReversed && pos == other.pos && refName == other.refName;
    }
};
// ReadFreq 结构体
struct ReadFreq {
    int freq;
    bam1_t* b;
    ReadFreq() : freq(0), b(nullptr) {}

    //int score;
    ReadFreq(bam1_t* b)  {
        this->freq=1;
        this->b=b;
    }
    void merge(bam1_t* b2){
        if (get_avg_qual(b2)> get_avg_qual(b)){
            this->b=b2;
        }
        freq++;
    }

};

// Hash 函数
namespace std {
    template <>
    struct hash<Alignment> {
        size_t operator()(const Alignment &a) const {
            return hash<string>()(a.refName) ^ hash<int>()(a.pos) ^ hash<bool>()(a.isReversed);
        }
    };
}
int main(int argc,char *argv[]){
    samFile *bam_in= sam_open(argv[1],"r");
    cout<<argv[1]<<endl;
    sam_hdr_t *bam_header= sam_hdr_read(bam_in);
    hts_set_threads(bam_in,20);
    bam1_t *b=bam_init1();
    double start_time,end_time;
    start_time=omp_get_wtime();
    char sep='_';
    int i=1;
    int unmapped=0;
    unordered_map<Alignment, unordered_map<string , ReadFreq>> align;
    while (sam_read1(bam_in,bam_header,b)>=0){
        //判断是否是无效印迹
        if (b->core.flag&BAM_FUNMAP){
            unmapped++;
            continue;
        }
        bool readNegativeFlag=b->core.flag&BAM_FREVERSE;
        string q_name= bam_get_qname(b);
        string umi= extractUMI(q_name,sep);
        Alignment alignment{readNegativeFlag,
                            (readNegativeFlag ? get_unclipped_end(b) : get_unclipped_start(b)),
                            get_reference_name(b,bam_header)};

        // 确保对齐位置已存在
        if (!align.count(alignment)) {
            align[alignment] = unordered_map<string, ReadFreq>();
        }

        // 处理 UMI 计数
        auto &umiMap = align[alignment];

        if (umiMap.count(umi)) {
            umiMap[umi].merge(b);
        } else {
            umiMap[umi] = ReadFreq(b);
        }
        i++;

    }
    cout<<"sum is "<<i<<",cost time:"<<omp_get_wtime()-start_time<<endl;
    cout<<"unmapped number:"<<unmapped<<endl;
    bam_destroy1(b);
    bam_hdr_destroy(bam_header);
    sam_close(bam_in);
//    hts_idx_destroy();
}

/**
 * 从read_name/q_name中提取umi
 * 已验证，结果相同
 * @param readName
 * @param sep umi分格符
 * @return
 */
std::string extractUMI(const std::string& readName, char sep) {
    // 找到分隔符的位置
    size_t sepPos = readName.find(sep);

    // 如果找不到分隔符，返回空字符串
    if (sepPos == std::string::npos) {
        return "";
    }

    // 提取分隔符后的部分，即 UMI
    return readName.substr(sepPos + 1);  // `sepPos + 1` 是跳过分隔符字符
}
/**
 * 获取平均质量分数，原版就是int，可以考虑换为float提高精度
 * 已验证，结果相同
 * @param b
 * @return
 */
int get_avg_qual(bam1_t* b){
    uint8_t *qual= bam_get_qual(b);
    int read_length=b->core.l_qseq;
    long sum_qual=0;
    for (int j = 0; j < read_length; ++j) {
        sum_qual+=qual[j];
    }
    int avg_qual = (sum_qual) / read_length;
    return avg_qual;
}
/**
 * 获取未裁剪的开始，用于确定位置
 * @param b
 * @return
 */
int get_unclipped_start(bam1_t *b) {
    int unclipped_start = b->core.pos+1;
    uint32_t *cigar = bam_get_cigar(b);
    int i;
    for (i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            unclipped_start -= len;
        } else {
            break;
        }
    }
    return unclipped_start;
}
/**
 * 获取未裁剪的末尾，用于确定位置
 * 已验证
 * @param b
 * @return
 */
int get_unclipped_end(bam1_t *b) {
    int unclipped_end = bam_endpos(b);
    uint32_t *cigar = bam_get_cigar(b);

    // 遍历 CIGAR 字符串，从后向前处理
    for (int i = b->core.n_cigar - 1; i >= 0; i--) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        // 软裁剪和硬裁剪时，调整结束位置
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            unclipped_end += len;
        } else {
            break;  // 其他操作，停止调整
        }
    }
    return unclipped_end;
}
/**
 * 获取引用名称，大部分为1，用于确定位置
 * 已验证
 * @param b
 * @param header
 * @return
 */
const char* get_reference_name(bam1_t *b, sam_hdr_t *header) {
    // 从 BAM 记录获取参考序列的索引（从 0 开始）
    int reference_index = b->core.tid;

    // 如果参考索引无效，返回 NULL
    if (reference_index < 0) {
        return NULL;
    }

    // 使用 BAM 文件头来获取参考名称
    return header->target_name[reference_index];
}