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
int main(int argc,char *argv[]){
    samFile *bam_in= sam_open(argv[1],"r");
    cout<<argv[1]<<endl;
    bam_hdr_t *bam_header= sam_hdr_read(bam_in);
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
//        uint8_t *umi_tag= bam_aux_get(aln,"RX");

        cout<<"read name:"<<extractUMI(bam_get_qname(b),sep)<<endl;
        if (i==100){
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