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
#include <vector>
#include <bitset>
using namespace std;

vector<bam1_t*>b_array(140000000);
using umi_type=bitset<64>  ;
int get_avg_qual(bam1_t* b);
// Initialize the static member
umi_type extractUMI(const std::string& readName, char sep);
int get_unclipped_start(bam1_t *b) ;
int get_unclipped_end(bam1_t *b);
//a,c,g,t之间的编码不同的位数都是2
const int ENCODING_DIST=2;
const std::unordered_map<char, uint8_t> ENCODING_MAP{
        {'A', 0b000},
        {'T', 0b101},
        {'C', 0b110},
        {'G', 0b011},
        {'a', 0b000},
        {'b', 0b101},
        {'c', 0b110},
        {'g', 0b011}
};

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
    unsigned int b;
    ReadFreq() : freq(0), b(0) {}

    //int score;
    ReadFreq(unsigned int b)  {
        this->freq=1;
        this->b=b;
    }
    void merge(unsigned int b2){
        if (get_avg_qual(b_array[b2])> get_avg_qual(b_array[this->b])){
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
    omp_set_num_threads(50);

    samFile *bam_in= sam_open(argv[1],"r");
    cout<<argv[1]<<endl;
    sam_hdr_t *bam_header= sam_hdr_read(bam_in);
    hts_set_threads(bam_in,50);

    double start_time,end_time;
    start_time=omp_get_wtime();
    char sep='_';
    unsigned int read_count=0;
    int unmapped=0;
    unordered_map<Alignment, unordered_map<umi_type , ReadFreq>> align;
//    vector<bam1_t*>b_array(100000000);
//    vector<bam1_t*>header_array(100000000);
//    bam1_t *b=bam_init1();
    while (true){
        bam1_t *b=bam_init1();
        if (sam_read1(bam_in,bam_header,b)<0){
            break;
        }
        b_array[read_count]=b;
//        header_array[i]=bam_header;
        //判断是否是无效印迹
        if (b->core.flag&BAM_FUNMAP){
            unmapped++;
            continue;
        }

        bool readNegativeFlag=b->core.flag&BAM_FREVERSE;
        string q_name= bam_get_qname(b);
        umi_type umi= extractUMI(q_name,sep);
        Alignment alignment{readNegativeFlag,
                            (readNegativeFlag ? get_unclipped_end(b) : get_unclipped_start(b)),
                            get_reference_name(b,bam_header)};


        // 确保对齐位置已存在
        if (!align.count(alignment)) {
            align[alignment] = unordered_map<umi_type, ReadFreq>();
        }
//        cout<<get_unclipped_start(b)<<endl;
        // 处理 UMI 计数
        auto &umiMap = align[alignment];

        auto [iter, inserted] = umiMap.emplace(umi, read_count);
        if (!inserted) {
            iter->second.merge(read_count);
        }


        read_count++;
    }
    cout<<"map is over"<<endl;
    cout<<"unmapped number:"<<unmapped<<endl;
    cout<<"map over,cost time:"<<omp_get_wtime()-start_time<<endl;
    unsigned int alignPosCount = align.size();
    unsigned int avgUMICount = 0,maxUMICount = 0,dedupedCount = 0;
    vector<std::pair<Alignment, unordered_map<umi_type , ReadFreq>>> entries(align.begin(),align.end());

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < entries.size(); ++i) {
        auto& [ali, umiMap]=entries[i];
        vector<int>deduped(umiMap.size());
        vector<pair<umi_type,ReadFreq>> freqs(umiMap.begin(),umiMap.end());
        sort(freqs.begin(), freqs.end(),
             [](const pair<umi_type, ReadFreq>& a, const pair<umi_type, ReadFreq>& b) {
                 return a.second.freq > b.second.freq;  // 降序排序
             }
        );

    }
//    bam_destroy1(b);
    cout<<"sum is "<<read_count<<",cost time:"<<omp_get_wtime()-start_time<<endl;
    #pragma omp parallel for
    for (int i = 0; i < read_count; ++i) {
        bam_destroy1(b_array[i]);
    }
    bam_hdr_destroy(bam_header);
    sam_close(bam_in);
//    free(b_array);
    cout<<"close and destroy,cost time:"<<omp_get_wtime()-start_time<<endl;
    return 0;
//    hts_idx_destroy();
}

/**
 * 从read_name/q_name中提取umi
 * 已验证，结果相同
 * @param readName
 * @param sep umi分格符
 * @return
 */
umi_type extractUMI(const std::string& readName, char sep) {
    umi_type result;
    // 找到分隔符的位置
    int sepPos = readName.find(sep);
    int bit_pos=0;
//    // 如果找不到分隔符，返回空字符串
    if (sepPos == std::string::npos) {
        return result;
    }
    string umi_str=readName.substr(sepPos + 1);  // `sepPos + 1` 是跳过分隔符字符
//    cout<<"umi: "<<umi_str<<endl;

    for (const char c : umi_str) {
        const auto it = ENCODING_MAP.find(c);
        if (it == ENCODING_MAP.end()) {
            return result;  // 遇到無效字符返回空bitset
        }

        const uint8_t code = it->second;
        // 將3bit編碼按高位到低位依次寫入bitset
        for (int shift = 2; shift >= 0; --shift) {
            if (bit_pos >= result.size()) {
                return result;  // 超過最大長度時截斷
            }
            result.set(bit_pos++, (code >> shift) & 0x1);
        }
    }
//    cout<<"umi bitset: "<<result.to_string()<<endl;

    // 提取分隔符后的部分，即 UMI
    return result;
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
int umi_dist(bitset<64> &a,bitset<64>& b){
    return (a^b).count()/2;
}
