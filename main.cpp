//
// Created by laoyao on 2025/1/27.
//
#include "htslib/sam.h"

// #include "sam.h"

#include <iostream>
#include "hisat_3n_table.h"
#include "ThreadPool.h"
#include <omp.h>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <algorithm>
// #include "thread_pool.h"
#include "htslib/thread_pool.h"
#include <unistd.h>
#include <fstream>
#include <list>
#include <atomic>
#include <queue>
#include <getopt.h>
#include <cmath>
using namespace std;

//using umi_type=bitset<64>  ;
using umi_type=string  ;
using score_type=double;
using adj_type=vector<umi_type>;
score_type get_avg_qual(bam1_t* b);
// Initialize the static member
umi_type extractUMI(const std::string& readName);
int64_t get_unclipped_start(bam1_t *b) ;
int64_t get_unclipped_end(bam1_t *b);
int unmapped=0;
unsigned int alignPosCount=0 ;
unsigned int avgUMICount=0 ,maxUMICount=0,dedupedCount=0 ;
//near的乘值
double percentage = 0.5f;//p
sam_hdr_t *bam_header;
//near的系数
int k=1; // k
bool read_done= false;
const unsigned int CHUNK_SIZE=1000000;
alignas(64) unsigned int total_read_count=0;
unsigned int chunk_N=0;
bam1_t *** b_array;
bool* consumer_flags;
int num_of_consumer=33;
int num_of_hts=34;
int thread_num = 16;
//是否在运行时回收内存
bool save_memory= false; //-save
char sep='_';
//a,c,g,t之间的编码不同的位数都是2
const int ENCODING_DIST=2;
string refFileName;

bool filter_rlen= false; //rlen


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
struct laoyaoAlignment {
    bool isReversed;
    int64_t pos;
    int32_t tid;
    // string refName;
    bool operator==(const laoyaoAlignment &other) const {
        return isReversed == other.isReversed && pos == other.pos && tid == other.tid;
    }
};
void freeBam1_t(unsigned index);
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
            if (save_memory) {
                freeBam1_t(this->b);
            }
            this->b=b2;
            this->score=b2_score;
        } else{
            if (save_memory) {
                freeBam1_t(b2);
            }
        }
        freq++;
    }

};
//用于向消费者传输数据
vector<pair<laoyaoAlignment , unsigned int >>** processor_datas;
vector<pair<laoyaoAlignment , unsigned int >>** consumer_local_datas;
// Hash 函数
namespace std {
    template <>
    struct hash<laoyaoAlignment> {
        size_t operator()(const laoyaoAlignment &a) const {
            return hash<int32_t>()(a.tid) ^ hash<int64_t>()(a.pos) ^ hash<bool>()(a.isReversed);
        }
    };
}
unsigned long getHash(const laoyaoAlignment& ali){
    return (hash<bool>()(!ali.isReversed)^hash<int32_t>()(ali.tid) ^ hash<int64_t>()(ali.pos)) ;
}
//unordered_map<Alignment, unordered_map<umi_type , ReadFreq*>> alignmentMap;
void visitAndRemove(const umi_type & u,unordered_map<umi_type, unordered_set<umi_type>>& adj,
                    unordered_set<umi_type>& visited);
void visitAndRemoveV2(const umi_type & u,
                      unordered_map<umi_type, adj_type>& adj, unordered_set<umi_type>& visited);
int umi_dist(const umi_type &a,const umi_type& b);
//unsigned int mapSectionChunk(sam_hdr_t *bam_header,unsigned int chunk_n,unsigned int read_size);
unsigned int produceChunk(sam_hdr_t *bam_header,unsigned int chunk_n,unsigned int read_size);
void consumeOne(const laoyaoAlignment& alignment,unsigned int pos,
                unordered_map<laoyaoAlignment, unordered_map<umi_type , ReadFreq*>>& local_align,
                unordered_map<laoyaoAlignment, unordered_map<umi_type,adj_type>>& align_adjs);
bam1_t* getBam1_t(unsigned int index){
    return b_array[index/CHUNK_SIZE][index%CHUNK_SIZE];
}



bool cmp(unsigned int i1,unsigned int i2){
//    return (getBam1_t(i1)->core.pos)<(getBam1_t(i2)->core.pos);
    const bam1_t *rec_a = getBam1_t(i1);
    const bam1_t *rec_b = getBam1_t(i2);
    // 比较染色体ID（tid）
    // if (rec_a->core.tid != rec_b->core.tid) {
    //     return rec_a->core.tid <rec_b->core.tid;
    // }

    // 比较起始位置（pos）
    if (rec_a->core.pos != rec_b->core.pos) {
        return rec_a->core.pos <rec_b->core.pos;
    }

    // 比对方向：正链（0x10标志位未设置）优先于负链（0x10标志位设置）
    int is_rev_a = bam_is_rev(rec_a);
    int is_rev_b = bam_is_rev(rec_b);
    return is_rev_a < is_rev_b;
}

bool cmp(const laoyaoAlignment rec_a,const laoyaoAlignment rec_b){
    //    return (getBam1_t(i1)->core.pos)<(getBam1_t(i2)->core.pos);


    // 比较起始位置（pos）
    if (rec_a.pos != rec_b.pos) {
        return rec_a.pos < rec_b.pos;
    }

    if (rec_a.isReversed != rec_b.isReversed) {
        return rec_a.isReversed < rec_b.isReversed;
    }
    return rec_a.tid < rec_b.tid;;
}
void freeBam1_t(unsigned index) {
//    free(b_array[index / CHUNK_SIZE][index % CHUNK_SIZE]->data);
//    free(b_array[index / CHUNK_SIZE][index % CHUNK_SIZE]);
    bam_destroy1(b_array[index / CHUNK_SIZE][index % CHUNK_SIZE]);
    b_array[index / CHUNK_SIZE][index % CHUNK_SIZE] = nullptr;
}
int32_t get_ref_length(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int32_t ref_len = 0;
    for (int i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op = cigar[i];
        char c = "MIDNSHP=XB"[op & BAM_CIGAR_MASK];
        int len = op >> BAM_CIGAR_SHIFT;
        if (c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X') {
            ref_len += len;
        }
    }
    return ref_len;
}
//bam1_t *bam_initV2(void)
//{
//    return (bam1_t*)malloc( sizeof(bam1_t));
//}

bool fileExist (string& filename) {
    ifstream file(filename);
    return file.good();
}

//参数解析
static const char *short_options = "i:o:t:p:k:fs"; 
static struct option long_options[] {
    {"input",  required_argument, 0, 'i'},
    {"output",  required_argument, 0, 'o'},
    {"thread", required_argument, 0, 't'},
    {"percentage", required_argument, 0, 'p'},
    {"k-near", required_argument, 0, 'k'},
    {"filter", no_argument, 0, 'f'},
    {"save-memory", no_argument, 0, 's'},
    {"ref", required_argument, 0, 'r' },
    {0, 0, 0, 0}
};



string inputFile, outputFilePath;

static void parseOption(int next_option, const char *optarg) {
    switch (next_option) {
        case 'i': {
            inputFile = optarg;
            if (!fileExist(inputFile)) {
                cout << "The input file is not exist" << endl;
                throw (1);
            }
            break;
        }
        case 'o': {
            outputFilePath = optarg;
            if(outputFilePath.empty()) {
                cout << "the outputFilePath is emtpy" << endl;
                throw(1);
            }
            break;
        }
        case 't': {
            try {
                thread_num = stoi(optarg);
            } catch(...) {
                cout << "the type of thread_num is integer" << endl;
                throw 1;
            }
            num_of_consumer = thread_num - 3;
            num_of_hts = thread_num - 1;
            break;
        }
        case 'p': {
            try {
                percentage = std::stof(optarg); 
            } catch(...) {
                cout << "the type of the percentage is float" << endl;
                throw 1;
            }
            break;
        }
        case 'k': {
            try {
                k = stoi(optarg);
            } catch (...) {
                cout << "the type of the k-near is integer" << endl;
            }
            break;
        }
        case 'f': {
            filter_rlen = true; //默认为false
            break;
        }
        case 's': {
            save_memory = true; //默认为false
            break;
        }
        case 'r': {
            refFileName = optarg;
            break;
        }

        default:
            cout << "please check the arguments" << endl;
            throw 1;
    } 
}
static void parseOptions(int argc, const char **argv) {
    int option_index = 0;
    int next_option;
    while(true) {
        next_option = getopt_long(argc, const_cast<char **>(argv), short_options, long_options, &option_index);
        if(next_option == -1)
            break;
        parseOption(next_option, optarg);
    }
    if(inputFile.empty()){
        cerr << "Error Non Input" << endl;
        exit(-1);
    }
    if(outputFilePath.empty()){
        cerr << "Error Non output" << endl;
        exit(-1);
    }

}

int main(int argc, const char** argv){
    try {
        parseOptions(argc, argv);
    } catch(...) {
        cerr << "some wrong happened when parse the options" << endl;
        return 1;
    }
    b_array=(bam1_t***) malloc(sizeof (bam1_t**)*100000);
    consumer_flags= (bool *)calloc(num_of_consumer, sizeof(bool));

    processor_datas= (vector<pair<laoyaoAlignment , unsigned int >>**)
            malloc(sizeof (vector<pair<laoyaoAlignment , unsigned int >>*)*num_of_consumer);

    consumer_local_datas=(vector<pair<laoyaoAlignment , unsigned int >>**)
            malloc(sizeof (vector<pair<laoyaoAlignment , unsigned int >>*)*num_of_consumer);
    omp_set_max_active_levels(2);
    omp_set_dynamic(0);

    for (int i = 0; i < num_of_consumer; ++i) {
//        vector<pair<Alignment , unsigned int >> p_v;
//        vector<pair<Alignment , unsigned int >> c_v;
        consumer_local_datas[i]=new vector<pair<laoyaoAlignment , unsigned int >>;
        processor_datas[i]=new vector<pair<laoyaoAlignment , unsigned int >>;
    }

    omp_set_num_threads(num_of_consumer);
    samFile *bam_in;
    bam_in = sam_open(inputFile.c_str(),"r");  //argument 1
    bam_header= sam_hdr_read(bam_in);
    htsThreadPool tpool = {NULL, 0};
    tpool.pool = hts_tpool_init(num_of_hts);
    if (tpool.pool) {
        hts_set_opt(bam_in, HTS_OPT_THREAD_POOL, &tpool);
    }
//    hts_set_threads(bam_in,20);
//    hts_set_threads(bam_out,20);
    double start_time;
    start_time=omp_get_wtime();
    avgUMICount = 0,maxUMICount = 0,dedupedCount = 0;

    unsigned int read_count=0;
    int init_done=0;
    unsigned int p_total=0;

    //按照染色体分好的数组
    auto final_write=(vector<unsigned int>***)malloc(sizeof(vector<unsigned int>**)*num_of_consumer);


#pragma omp parallel sections num_threads(3)
    {
        //初始化读入内存
        #pragma omp section
        {
            unsigned int p_read=0;
            auto ** b_array_local= (bam1_t** )malloc(CHUNK_SIZE*sizeof(bam1_t*));
            unsigned int p_chunk_n=0;
            while (true) {
                bam1_t* b= bam_init1();
                if (sam_read1(bam_in, bam_header, b) < 0) {
                    //执行完成
                    p_total+=p_read;
                    if (p_read>0){
                        b_array[p_chunk_n]=b_array_local;
                    }
                    total_read_count=p_total;
                    #pragma omp atomic seq_cst
                    init_done++;
                    if (p_read>0) {
                        #pragma omp atomic seq_cst//acquire
                        chunk_N++;
                    }
                    #pragma omp atomic seq_cst
                    init_done++;
                    bam_destroy1(b);
                    cout << "init vector over,cost time:" << omp_get_wtime() - start_time << endl;
                    cout<<total_read_count<<endl;
                    break;
                }
                b_array_local[p_read]=b;
                p_read++;
                //每一定间隔发送数据
                if (p_read%CHUNK_SIZE==0) {
                    //将这块数据放入
                    b_array[p_chunk_n++]=b_array_local;
                    p_total+=CHUNK_SIZE;
                    #pragma omp atomic seq_cst
                    chunk_N++;
                    p_read=0;
                    b_array_local=(bam1_t** )malloc(sizeof(bam1_t*)*CHUNK_SIZE);
//                    cout<<"send chunk"<<endl;
                }
            }
        }
        //
        #pragma omp section
        {

            unsigned int current_n=0;
            unsigned int toN=0;
            cout<<"producer begin "<<endl;
        //进入第一步map阶段
            while (true) {
                #pragma omp atomic read seq_cst
                toN=chunk_N;
                if (init_done) {
                    while (init_done!=2){
                        usleep(100);
                    }
                    #pragma omp atomic read seq_cst
                    toN=chunk_N;
                    for (unsigned int i = current_n; i < toN-1; ++i) {
                        read_count+= produceChunk(bam_header,i,CHUNK_SIZE);
                    }
                    read_count+=produceChunk(bam_header,toN-1,(total_read_count-1)%CHUNK_SIZE+1);
                    //确保生产者的所有数据被发送给消费者
                    while (true) {
                        int number_rest=num_of_consumer;
                        for (int i = 0; i < num_of_consumer; ++i) {
                            if(consumer_local_datas[i]== nullptr||consumer_local_datas[i]->empty()){
                                number_rest--;
                                continue;
                            }
                            //内部无数据
                            if (!consumer_flags[i]) {
                                processor_datas[i] = consumer_local_datas[i];
                                #pragma omp atomic write seq_cst
                                //成功放入后
                                consumer_flags[i] = true;
                                number_rest--;
                                consumer_local_datas[i] = nullptr;
                            }
                        }
                        if (number_rest==0){
                            break;
                        }
                    }
                    cout<<"all producer local data clear"<<endl;
                    break;
                }
                if (toN==current_n){
                    usleep(1000);
                    continue;
                }
//                cout<<"receive chunk"<<toN<<endl;

                for (unsigned int i = current_n; i < toN; ++i) {
                    read_count+= produceChunk(bam_header,i, CHUNK_SIZE);
                }
                current_n=toN;
                //不再有新的，处理完剩余的退出
            }



            #pragma omp atomic write seq_cst
            read_done= true;
        }
        #pragma omp section
        {
            #pragma omp parallel num_threads(num_of_consumer) reduction(+:dedupedCount,avgUMICount,alignPosCount) reduction(max:maxUMICount)
            {
                int id = omp_get_thread_num();
                //将每个alignment 下面的所有位置索引存储起来,并提前分组
                auto* local_sort_Alignment=(vector<unsigned int>**)malloc(sizeof(vector<unsigned int>**)*25);
                for (int i = 0; i < 25; ++i) {
                    local_sort_Alignment[i] = new vector<unsigned int>;
                }

//                cout<<id<<" begin"<<endl;
                unordered_map<laoyaoAlignment, unordered_map<umi_type, ReadFreq *>> local_align;

                unordered_map<laoyaoAlignment, unordered_map<umi_type, adj_type>> align_adjs;
                vector<unsigned int> poses;
                vector<pair<laoyaoAlignment , unsigned int >> process_data;
                //获取信号，表示有新的
                while (true) {
                    if (consumer_flags[id]) {
                        process_data=*processor_datas[id];
                        #pragma atomic write seq_cst
                        consumer_flags[id]= false;
                        for (const auto &[ali, pos]: process_data) {
                            consumeOne(ali, pos, local_align, align_adjs);
                        }
                        //等待一段事件再次获取
                    } else {
                        usleep(100);
                    }

                    if (read_done) {
                        //再一次执行，确保处理完全
                        if (consumer_flags[id]) {
                            process_data=*processor_datas[id];
                            for (const auto &[ali, pos]: process_data) {
                                consumeOne(ali, pos, local_align, align_adjs);
                            }
                        }
                        alignPosCount=local_align.size();
                        //处理完成退出线程
                        break;
                    }
                }

//                vector<unsigned int> dedupeds;
                //处理完成进入后处理阶段
                for (auto &tackle_map: align_adjs) {

                    unordered_map<umi_type, ReadFreq *> umi_read = local_align[tackle_map.first];

                    vector<pair<umi_type, ReadFreq *>> freqs(umi_read.begin(), umi_read.end());
                    //本地当前的alignemnt
                    laoyaoAlignment current_aligment=tackle_map.first;

//
                    sort(freqs.begin(), freqs.end(),
                         [](const pair<umi_type, ReadFreq *> &a, const pair<umi_type, ReadFreq *> &b) {
                             return a.second->freq > b.second->freq;  // 降序排序
                         }
                    );

                    unordered_map<umi_type, adj_type> umiMap = tackle_map.second;
                    //再次删除,根据freq删除
                    for (auto &umi_adjs: umiMap) {
                        int maxFreq = static_cast<int>(((umi_read[umi_adjs.first]->freq) + 1) * percentage);
                        adj_type &to_remove = umi_adjs.second;
                        auto new_end = std::remove_if(
                                to_remove.begin(),
                                to_remove.end(),
                                [&](const auto &elem) {
                                    return (umi_read[elem]->freq) > maxFreq;
                                }
                        );
                        to_remove.erase(new_end, to_remove.end());
                    }
                    unordered_set<umi_type> visited;

                    //存储这个alignment下的所有结果
                    auto the_rest_of_alignment=new vector<unsigned int>;
                    for (const auto &freq: freqs) {

                        if (visited.count(freq.first) == 0) {
                            visitAndRemoveV2(freq.first, umiMap, visited);
                            //放入bed
                            the_rest_of_alignment->emplace_back(freq.second->b);
                            // local_sort_chunk[get_bucket_index(freq.second->b)]->push_back(freq.second->b);
                            dedupedCount++;
                        }
                    }

                    if (current_aligment.tid>23) {
                        local_sort_Alignment[24]->insert(local_sort_Alignment[24]->end(),the_rest_of_alignment->begin(), the_rest_of_alignment->end());
                    }else {
                        local_sort_Alignment[current_aligment.tid]->insert(local_sort_Alignment[current_aligment.tid]->end(),the_rest_of_alignment->begin(), the_rest_of_alignment->end());
                    }

                    avgUMICount += umi_read.size();
                    maxUMICount = std::max(maxUMICount, (unsigned int) umi_read.size());

                    //开始写入
                }
                //将本地分组结果写入
                final_write[id]=local_sort_Alignment;

            }
        }
    }
    cout<<"total read count"<<total_read_count<<endl;

    cout<<"unmapped number:"<<unmapped<<endl;
    cout<<"map is over"<<" current time:"<<omp_get_wtime()-start_time<<endl;



    auto merged_positions=(vector<unsigned int>**)malloc(sizeof(vector<unsigned int>*) * 25);;

    // 预计算总元素数以优化内存分配
    #pragma omp parallel num_threads(25)
    {

        int tid = omp_get_thread_num();
        if (tid == 0) {
            cout<<" tid 0000000000000000"<<endl;
        }
        auto current_pos=new vector<unsigned int>;
        for (int i = 0; i < num_of_consumer; ++i) {
            current_pos->insert(current_pos->end(),final_write[i][tid]->begin(), final_write[i][tid]->end());
        }
        if (tid<24) {
            sort(current_pos->begin(), current_pos->end(),
            [](unsigned int i1, unsigned int i2) {
            // 内联优化版比较逻辑
                const bam1_t* a = getBam1_t(i1);
                const bam1_t* b = getBam1_t(i2);
            // 优化分支预测顺序
                return (a->core.pos != b->core.pos) ? (a->core.pos < b->core.pos) :
                    (bam_is_rev(a) < bam_is_rev(b));
                });
        }else {
            sort(current_pos->begin(), current_pos->end(),
            [](unsigned int i1, unsigned int i2) {
                // 内联优化版比较逻辑
                    const bam1_t* a = getBam1_t(i1);
                    const bam1_t* b = getBam1_t(i2);
                // 优化分支预测顺序
                    return (a->core.tid != b->core.tid) ? (a->core.tid < b->core.tid) :
                   (a->core.pos != b->core.pos) ? (a->core.pos < b->core.pos) :
                   (bam_is_rev(a) < bam_is_rev(b));
            });
        }
        merged_positions[tid]=current_pos;

    }




    ThreadPool h3tThreadPool(60);

    std::string base_output_path = outputFilePath; //argument 2

    string outputPath=base_output_path;
    // cout<<"merge done start to write:"<<offset<<" current time:"<<omp_get_wtime()-start_time<<endl;
    //对0号
    int case_name=0;

    outputPath = base_output_path;

    for (int i = 0; i <25 ; ++i) {

        //必须在这里提交任务
        h3tThreadPool.enqueue([merged_positions, i, outputPath, case_name](){
            hisat_3n_table(merged_positions[i]->data(), 0, merged_positions[i]->size()-1, outputPath + "unfiltered_unique." + "part" + to_string(case_name+1) + ".tsv", true, false, false);
        });
        h3tThreadPool.enqueue([merged_positions, i, outputPath, case_name](){
            hisat_3n_table(merged_positions[i]->data(), 0, merged_positions[i]->size()-1, outputPath + "unfiltered_multi." + "part" + to_string(case_name+1) + ".tsv", false, true, false);
        });
        h3tThreadPool.enqueue([merged_positions, i, outputPath, case_name](){
            hisat_3n_table(merged_positions[i]->data(), 0, merged_positions[i]->size()-1, outputPath + "filtered_unique." + "part" + to_string(case_name+1) + ".tsv", true, false, true);
        });
        h3tThreadPool.enqueue([merged_positions, i, outputPath, case_name](){
            hisat_3n_table(merged_positions[i]->data(), 0, merged_positions[i]->size()-1, outputPath + "filtered_multi." + "part" + to_string(case_name+1) +  ".tsv", false, true, true);
        });
        case_name++;


    }

    h3tThreadPool.stop_pool();


    cout<<"sum is "<<read_count<<",cost time:"<<omp_get_wtime()-start_time<<endl;
    cout<<"Number of unremoved reads\t" <<read_count<<endl;
    cout<<"Number of unique alignment positions\t" <<alignPosCount<<endl;
    cout<<"Average number of UMIs per alignment position\t" <<((double)avgUMICount / alignPosCount)<<endl;
    cout<<"Max number of UMIs over all alignment positions\t" << maxUMICount<<endl;
    cout<<"Number of reads after deduplicating\t" << dedupedCount<<endl;
    sam_close(bam_in);

    bam_hdr_destroy(bam_header);
//    if (tpool.pool) {
//        hts_tpool_destroy(tpool.pool);  // Must be done AFTER sam_close()
//    }
    cout<<"close and destroy,cost time:"<<omp_get_wtime()-start_time<<endl;
    exit(0);
//    bam_destroy1(b);
    #pragma omp parallel for schedule(dynamic) num_threads(num_of_hts)
    for (int i = 0; i < chunk_N; ++i) {
        unsigned int last=CHUNK_SIZE;
        if (i==CHUNK_SIZE-1){
            last=total_read_count%CHUNK_SIZE;
        }
        bam1_t ** to_destroy=b_array[i];
        for (int j = 0; j < last; ++j) {
            if (!save_memory) {
                bam_destroy1(to_destroy[j]);
            } else {
                if (to_destroy[i] != nullptr) {
                    bam_destroy1(to_destroy[j]);
                }
            }
        }
        free(to_destroy);
    }
//    free(b_array);
    return 0;
//    hts_idx_destroy();
}
/**
 * 根据alignment分发数据
 * @param bam_header 头
 * @param chunk_n 处理第几快
 * @param read_size 处理多少个，一般为CHUNK
 * @return
 */
unsigned int produceChunk(sam_hdr_t *bam_header,unsigned int chunk_n,unsigned int read_size){
    int local_read_counts=0;
    bam1_t** chunk=b_array[chunk_n];
    for (int p_read_count = 0; p_read_count <read_size; ++p_read_count) {
        bam1_t *b = chunk[p_read_count];
        //判断是否是无效印迹s
        if (b->core.flag & BAM_FUNMAP) {
            unmapped++;
            continue;
        }
        bool readNegativeFlag = b->core.flag & BAM_FREVERSE;
        //获取这条记录的位置
        laoyaoAlignment alignment{readNegativeFlag,
                            (readNegativeFlag ? get_unclipped_end(b) : get_unclipped_start(b)),
                            b->core.tid};
        unsigned long pos_hash= getHash(alignment)%num_of_consumer;
        consumer_local_datas[pos_hash]->emplace_back(alignment,chunk_n*CHUNK_SIZE+p_read_count);
        local_read_counts++;

    }
    for (int i = 0; i < num_of_consumer; ++i) {

        //内部还有数据，
        if (!consumer_flags[i]){
//            cout<<"ready to put "<<i<<endl;
            processor_datas[i]=consumer_local_datas[i];
//            cout<<"put sucess"<<i<<endl;
            #pragma omp atomic write seq_cst
            //成功放入后
            consumer_flags[i]= true;
//            vector<pair<Alignment,unsigned int >> vec_tmp;
            consumer_local_datas[i]=new vector<pair<laoyaoAlignment,unsigned int >>;
        }

    }
    return local_read_counts;
}


void consumeOne(const laoyaoAlignment& alignment,unsigned int pos,
                unordered_map<laoyaoAlignment, unordered_map<umi_type , ReadFreq*>>& local_align,
unordered_map<laoyaoAlignment, unordered_map<umi_type,adj_type>>& align_adjs){
    bam1_t * b= getBam1_t(pos);
    if (filter_rlen){
        if (get_ref_length(b)>=100000){
            if (save_memory){
                freeBam1_t(pos);
            }
            return;
        }
    }
    string q_name= bam_get_qname(b);
    umi_type umi= extractUMI(q_name);
    score_type score= get_avg_qual(b);

    auto &umiMap = local_align[alignment];
    auto it_read = umiMap.find(umi);
    if (it_read != umiMap.end()) {
        //原来就已经存在，保存质量最好的一条
        it_read->second->merge(pos,score );
    } else {
        //原来不存在合并
        //将结果放入对应adj,根据dist筛选
        unordered_map<umi_type,adj_type>adjs=align_adjs[alignment];
        adj_type one_adj;
        //提前计算每个umi之间的距离
        for (auto & iter : umiMap) {
            int dist=umi_dist(iter.first,umi);
            if (dist<=k){
                adjs[iter.first].push_back(umi);
                one_adj.push_back(iter.first);
            }
        }
        umiMap.insert({umi, new ReadFreq(pos, score)});
        adjs[umi]=one_adj;
        align_adjs[alignment]=adjs;
    }

}

/**
 * 从read_name/q_name中提取umi
 * 已验证，结果相同
 * @param readName
 * @param sep umi分格符
 * @return
 */
umi_type extractUMI(const std::string& readName) {
    umi_type result;
    // 找到分隔符的位置
    unsigned long sepPos = readName.find(sep);
//    // 如果找不到分隔符，返回空字符串
    if (sepPos == std::string::npos) {
        return result;
    }
    result=readName.substr(sepPos + 1);

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
    score_type avg_qual = static_cast<score_type>(sum_qual) / read_length;
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
int umi_dist(const umi_type &a,const umi_type& b){
    int result=0;
//    return (a^b).count()/ENCODING_DIST;
    for (int i = 0; i < a.size(); ++i) {
        if (a.at(i)!=b.at(i)){
            result++;
        }
    }
    return result;
}
unordered_set<umi_type> near(const vector<pair<umi_type,ReadFreq*>>& freqs,const umi_type& umi, int maxFreq){
    unordered_set<umi_type>res;
    for (auto & freq : freqs) {
        umi_type o=freq.first;
        int f=freq.second->freq;
        int dist= umi_dist(umi,o);
        if (dist<=k&&(f<=maxFreq)){
            if (res.count(o)!=0){
                cout<<" repeat!!!!!!!!!!!!"<<endl;
            }
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

void visitAndRemoveV2(const umi_type & u,
                     unordered_map<umi_type, adj_type>& adj, unordered_set<umi_type>& visited) {
    if (visited.count(u)!=0) return;
    // 获取邻接的节点
    const auto& neighbors = adj[u];
    visited.insert(u);
    for (const auto& v : neighbors) {
        if (u == v) continue;
        visitAndRemoveV2(v, adj, visited);
    }
}