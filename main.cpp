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
using namespace std;

// Initialize the static member
std::string extractUMI(const std::string& readName, char sep);
int get_unclipped_start(bam1_t *b) ;
int get_unclipped_end(bam1_t *b);

const char* get_reference_name(bam1_t *b, sam_hdr_t *header);

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
    while (sam_read1(bam_in,bam_header,b)>=0){
        //判断是否是无效印迹
        if (b->core.flag&BAM_FUNMAP){
            unmapped++;
            continue;
        }
        bool readNegativeFlag=b->core.flag&BAM_FREVERSE;
//        if (!readNegativeFlag) {
//            cout << "negative flag " << readNegativeFlag << endl;
//            cout<<"clip pos"<<get_unclipped_start(b)<<endl;
//            cout<<i<<endl;
////        std::cout << "length"<<read_length<< " average quality: " << avg_qual << std::endl;
////            cout << "read name: " << extractUMI(bam_get_qname(b), sep) << endl;
//        }
//        cout << "negative flag " << readNegativeFlag << endl;
//        if (readNegativeFlag){
//            cout << "clip pos; " << get_unclipped_end(b) << endl;
//        }else{
//            cout << "clip pos: " << get_unclipped_start(b) << endl;
//        }
//        cout<<"ref name"<<get_reference_name(b,bam_header)<<endl;
        if (strcmp(get_reference_name(b,bam_header),"1")!=0){
            cout << "negative flag " << readNegativeFlag << endl;
            if (readNegativeFlag){
                cout << "clip pos; " << get_unclipped_end(b) << endl;
            }else{
                cout << "clip pos: " << get_unclipped_start(b) << endl;
            }
            cout<<"ref name"<<get_reference_name(b,bam_header)<<endl;
            cout<<i<<endl;
            break;
        }

        if (i==50000000){
            break;
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
 * 已验证，与原版相比大小-1
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