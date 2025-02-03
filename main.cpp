//
// Created by laoyao on 2025/1/27.
//
#include "htslib/sam.h"
#include <iostream>
#include <omp.h>
#include <regex>
#include <bitset>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <algorithm>
//#include <cstring>
#include <algorithm>
#include "htslib/thread_pool.h"
#include <unistd.h>
using namespace std;

//using umi_type=bitset<64>  ;
using umi_type=string  ;
using score_type=int;
score_type get_avg_qual(bam1_t* b);
// Initialize the static member
umi_type extractUMI(const std::string& readName, char sep);
int64_t get_unclipped_start(bam1_t *b) ;
int64_t get_unclipped_end(bam1_t *b);
constexpr size_t POOL_SIZE = 135000000;
int unmapped=0;
const unsigned int CHUNK_SIZE=1000000;
alignas(64) unsigned int total_read_count=0;
bam1_t *** b_array;
char sep='_';
//a,c,g,t之间的编码不同的位数都是2
const int ENCODING_DIST=2;
const std::unordered_map<char, uint8_t> ENCODING_MAP{
        {'A', 0b000},
        {'T', 0b101},
        {'C', 0b110},
        {'G', 0b011},
        {'N', 0b100},
};

const char* get_reference_name(bam1_t *b, sam_hdr_t *header);
/**
 * 表示位置的结构体
 */
struct Alignment {
    bool isReversed;
    int64_t pos;
    string refName;
    bool operator==(const Alignment &other) const {
        return isReversed == other.isReversed && pos == other.pos && refName == other.refName;
    }
};
// ReadFreq 结构体
struct ReadFreq {
    int freq;
    unsigned int b;
    score_type score=0;
    ReadFreq() : freq(0), b(0) {}

    //int score;
    ReadFreq(unsigned int b,score_type score)  {
        this->freq=1;
        this->b=b;
        this->score=score;
    }
    void merge(unsigned int b2,score_type b2_score){
        if (b2_score > this->score){
            this->b=b2;
            this->score=b2_score;
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
unordered_map<Alignment, unordered_map<umi_type , ReadFreq*>> alignmentMap;
unordered_set<umi_type> near(const vector<pair<umi_type,ReadFreq*>>& freqs,umi_type umi,int k, int maxFreq);
void visitAndRemove(const umi_type & u,unordered_map<umi_type, unordered_set<umi_type>>& adj,
                    unordered_set<umi_type>& visited);
unsigned long umi_dist(umi_type &a,umi_type& b);
unsigned int mapSection(sam_hdr_t *bam_header,unsigned int start,unsigned int end);
bam1_t* getBam1_t(unsigned index){
    return b_array[index/CHUNK_SIZE-1][index%CHUNK_SIZE];
}

//bam1_t *bam_initV2(void)
//{
//    return (bam1_t*)malloc( sizeof(bam1_t));
//}

int main(int argc,char *argv[]){
    b_array=(bam1_t***) malloc(sizeof (bam1_t**)*CHUNK_SIZE);
    omp_set_num_threads(30);
    double percentage = 0.5f;
//    samFile *bam_in= sam_open("/root/input/SRR23538291.mRNA.genome.mapped.sorted.bam","r");
//    htsFile *bam_out= hts_open("/root/input/outputSRR23538290.mRNA.genome.mapped.sorted.dedup.bam","wb");
    samFile *bam_in= sam_open(argv[1],"r");
    htsFile *bam_out= hts_open(argv[2],"wb");
    sam_hdr_t *bam_header= sam_hdr_read(bam_in);
    if (sam_hdr_write(bam_out, bam_header) < 0) {
        cerr << "Error writing output." << endl;
        exit(-1);
    }
    cout<<argv[1]<<endl;
    htsThreadPool tpool = {NULL, 0};
    tpool.pool = hts_tpool_init(50);
    if (tpool.pool) {
        hts_set_opt(bam_in, HTS_OPT_THREAD_POOL, &tpool);
        hts_set_opt(bam_out, HTS_OPT_THREAD_POOL, &tpool);
    }
//    hts_set_threads(bam_in,20);
//    hts_set_threads(bam_out,20);
    double start_time,end_time;
    start_time=omp_get_wtime();
    int k=1;

    unsigned int read_count=0;
    alignas(64) bool init_done=false;
#pragma omp parallel sections num_threads(2)
    {
        //初始化读入内存
        #pragma omp section
        {
            cout << "bam init over,cost time:" << omp_get_wtime() - start_time << endl;
            unsigned int p_read=0;
            while (true) {
                b_array[total_read_count+p_read] = bam_init1();
                if (sam_read1(bam_in, bam_header, b_array[total_read_count+p_read]) < 0) {
                    //执行完成
                    #pragma omp atomic
                    init_done+=1;
                    #pragma omp atomic
                    total_read_count+=p_read;
                    bam_destroy1(b_array[total_read_count]);
                    cout << "init vector over,cost time:" << omp_get_wtime() - start_time << endl;
                    cout<<total_read_count<<endl;
                    break;
                }
                p_read++;
                //每一定间隔发送数据
                if (p_read%500000==0) {
                    #pragma omp atomic
                    total_read_count+=p_read;
                    p_read=0;
                }
            }
        }
        //
        #pragma omp section
        {
            usleep(2000000);
            unsigned int current_read=0;
            unsigned int toread=0;
        //进入第一步map阶段
            while (true) {
                #pragma omp atomic read
                toread=total_read_count;
                read_count += mapSection(bam_header,current_read,toread);
                current_read=toread;
                //不再有新的，处理完剩余的退出
                if (init_done) {
                    read_count+= mapSection(bam_header,current_read,total_read_count);
                    break;
                }
            }
        }
    }
//    exit(0);
//    #pragma omp parallel for reduction(+:read_count) reduction(+:unmapped)
    cout<<"total read count"<<total_read_count<<endl;

    cout<<"map is over"<<endl;
    cout<<"unmapped number:"<<unmapped<<endl;
    cout<<"map over,cost time:"<<omp_get_wtime()-start_time<<endl;

    unsigned int alignPosCount = alignmentMap.size();
    unsigned int avgUMICount = 0,maxUMICount = 0,dedupedCount = 0;
    vector<std::pair<Alignment, unordered_map<umi_type , ReadFreq*>>> entries(alignmentMap.begin(), alignmentMap.end());

    omp_set_num_threads(40);

#pragma omp parallel for schedule(dynamic) reduction(+:dedupedCount,avgUMICount) reduction(max:maxUMICount)
    for (unsigned int i = 0; i < entries.size(); ++i) {
        auto& [ali, umiMap]=entries[i];
        vector<ReadFreq*>deduped;
        //读取并按频率排序
        vector<pair<umi_type,ReadFreq*>> freqs(umiMap.begin(),umiMap.end());
//
        sort(freqs.begin(), freqs.end(),
             [](const pair<umi_type, ReadFreq*>& a, const pair<umi_type, ReadFreq*>& b) {
                 return a.second->freq > b.second->freq;  // 降序排序
             }
        );

        //构建邻链接边
        vector<unordered_set<umi_type>>adjIdx(freqs.size());
        unordered_map<umi_type,unordered_set<umi_type>>adj;
//        int near_number=0;
        for (int j = 0; j < freqs.size(); ++j) {
            adj[freqs[j].first]= near(freqs,freqs[j].first,k,
                            static_cast<int>((freqs[j].second->freq+1)*percentage));
//            near_number+=adj[freqs[j].first].size();
        }


        unordered_set<umi_type>visited;

        for(const auto& freq : freqs){
            if (visited.count(freq.first)==0){
                visitAndRemove(freq.first,adj,visited);
                deduped.push_back(freq.second);
            }
        }
        avgUMICount += umiMap.size();
        maxUMICount = std::max(maxUMICount, (unsigned int)umiMap.size());
        dedupedCount += deduped.size();
        //处理完成
        #pragma omp critical
        {
            for(ReadFreq* readFreq:deduped){
                if (sam_write1(bam_out,bam_header,b_array[readFreq->b])<0){
                    cerr << "Error writing output." << endl;
                    exit(-1);
                }
            }
        }
    }
//    bam_destroy1(b);
    cout<<"sum is "<<read_count<<",cost time:"<<omp_get_wtime()-start_time<<endl;
    cout<<"Number of unremoved reads\t" <<read_count<<endl;
    cout<<"Number of unique alignment positions\t" <<alignPosCount<<endl;
    cout<<"Average number of UMIs per alignment position\t" <<((double)avgUMICount / alignPosCount)<<endl;
    cout<<"Max number of UMIs over all alignment positions\t" << maxUMICount<<endl;
    cout<<"Number of reads after deduplicating\t" << dedupedCount<<endl;
    #pragma omp parallel for
    for (int i = 0; i < total_read_count; ++i) {
        bam_destroy1(b_array[i]);
    }
    bam_hdr_destroy(bam_header);
    sam_close(bam_in);
    if (tpool.pool) {
        hts_tpool_destroy(tpool.pool);  // Must be done AFTER sam_close()
    }
//    free(b_array);
    cout<<"close and destroy,cost time:"<<omp_get_wtime()-start_time<<endl;
    return 0;
//    hts_idx_destroy();
}
/**
 * map阶段，将所有位置，umi相同的记录聚合到一起
 * @param bam_header
 * @return
 */
unsigned int mapSection(sam_hdr_t *bam_header,unsigned int start,unsigned int end){
    int read_count=0;
    for (int p_read_count = start; p_read_count < end; ++p_read_count) {
        bam1_t* b=b_array[p_read_count];
        //判断是否是无效印迹s
        if (b->core.flag&BAM_FUNMAP){
            unmapped++;
            continue;
        }
        bool readNegativeFlag=b->core.flag&BAM_FREVERSE;
        string q_name= bam_get_qname(b);
        umi_type umi= extractUMI(q_name,sep);
//        cout<<q_name<<endl;
        //获取这条记录的位置
        Alignment alignment{readNegativeFlag,
                            (readNegativeFlag ? get_unclipped_end(b) : get_unclipped_start(b)),
                            get_reference_name(b,bam_header)};
        //获取这条记录的质量分数
        score_type score= get_avg_qual(b);
        #pragma omp critical
        {
            // 确保对齐位置已存在
            if (!alignmentMap.count(alignment)) {
                alignmentMap[alignment] = unordered_map<umi_type, ReadFreq *>();
            }
            // 处理 UMI 计数
            auto &umiMap = alignmentMap[alignment];
            auto it_read = umiMap.find(umi);
            if (it_read != umiMap.end()) {
                it_read->second->merge(p_read_count, score);
            } else {
                umiMap.insert({umi, new ReadFreq(p_read_count, score)});
            }
        }
        read_count++;
    }
    return read_count;
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
    result=readName.substr(sepPos + 1);
//    string umi_str=readName.substr(sepPos + 1);  // `sepPos + 1` 是跳过分隔符字符
//    cout<<"umi: "<<umi_str<<endl;

//    for (const char c : umi_str) {
//        const auto it = ENCODING_MAP.find(c);
//        if (it == ENCODING_MAP.end()) {
//            cout<<" un use it:"<<c<<endl;
//            return result;  // 遇到無效字符返回空bitset
//        }
//
//
//        const uint8_t code = it->second;
//        // 將3bit編碼按高位到低位依次寫入bitset
//        for (int shift = 2; shift >= 0; --shift) {
//            if (bit_pos >= result.size()) {
//                cout<<"being cuted "<<endl;
//                return result;  // 超過最大長度時截斷
//            }
//            result.set(bit_pos++, (code >> shift) & 0x1);
//        }
////        for (int shift = 0; shift < 3; ++shift) { // 修改循环方向
////            if (bit_pos >= result.size()) return result;
////            result.set(bit_pos++, (code >> shift) & 0x1);
////        }
//    }
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
score_type get_avg_qual(bam1_t* b){
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
int64_t get_unclipped_start(bam1_t *b) {
    int64_t unclipped_start = b->core.pos+1;
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
int64_t get_unclipped_end(bam1_t *b) {
    int64_t unclipped_end = bam_endpos(b);
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
        return nullptr;
    }

    // 使用 BAM 文件头来获取参考名称
    return header->target_name[reference_index];
}
/**
 * 计算两个umi中有多少个基因不相同
 */
unsigned long umi_dist(umi_type &a,umi_type& b){
    unsigned long result=0;
//    return (a^b).count()/ENCODING_DIST;
    for (int i = 0; i < a.size(); ++i) {
        if (a.at(i)!=b.at(i)){
            result++;
        }
    }
    return result;
}
unordered_set<umi_type> near(const vector<pair<umi_type,ReadFreq*>>& freqs,umi_type umi,int k, int maxFreq){
    unordered_set<umi_type>res;
    for (auto & freq : freqs) {
        umi_type o=freq.first;
        int f=freq.second->freq;
        unsigned long dist= umi_dist(umi,o);
        if (dist<=k&&(dist==0||f<=maxFreq)){
            res.insert(o);
        }
    }
    return res;

}
void visitAndRemove(const umi_type & u,
                    unordered_map<umi_type, unordered_set<umi_type>>& adj, unordered_set<umi_type>& visited) {
    if (visited.count(u)!=0) return;

    // 获取邻接的节点
    const auto& neighbors = adj[u];
    visited.insert(u);

    for (const auto& v : neighbors) {
        if (u == v) continue;
        visitAndRemove(v, adj, visited);
    }
}
