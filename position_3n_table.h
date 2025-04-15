#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>
#include "alignment_3n_table.h"
/**
 * store unique information for one base information with readID, and the quality.
 */
class uniqueID
{
public:
    unsigned long long readNameID;
    bool isConverted;
    char quality;
    bool removed;

    uniqueID(unsigned long long InReadNameID,
             bool InIsConverted,
             char &InQual)
    {
        readNameID = InReadNameID;
        isConverted = InIsConverted;
        quality = InQual;
        removed = false;
    }
};

/**
 * basic class to store reference position information
 */
class Position
{
public:
    string Pchromosome;          // reference chromosome name
    long long int Plocation;     // 1-based position
    char strand;                 // +(REF) or -(REF-RC)
    string convertedQualities;   // 每个字符表示该位置上已转换碱基的比对质量。
    string unconvertedQualities; // 每个字符表示该位置上未转换碱基的比对质量。
    vector<uniqueID> uniqueIDs;  // 每个值表示一个贡献了碱基信息的 readName。
                                 // readNameIDs 用于确保同一个 read 不会在同一个位置上贡献两次。

    void initialize()
    {
        Pchromosome.clear();
        Plocation = -1;
        strand = '?';
        convertedQualities.clear();
        unconvertedQualities.clear();
        vector<uniqueID>().swap(uniqueIDs);
    }

    Position()
    {
        initialize();
    };

    /**
     * return true if there is mapping information in this reference position.
     */
    bool empty()
    {
        return convertedQualities.empty() && unconvertedQualities.empty();
    }

    /**
     * set the chromosome, location (position), and strand information.
     */

    void set(string &inputChr, long long int inputLoc)
    {
        Pchromosome = inputChr;
        Plocation = inputLoc + 1;
    }

    void set(char inputStrand)
    {
        strand = inputStrand;
    }

    /**
     * binary search of readNameID in readNameIDs.
     * always return a index.
     * if cannot find, return the index which has bigger value than input readNameID.
     */
    int searchReadNameID(unsigned long long &readNameID, int start, int end)
    {
        if (uniqueIDs.empty())
        {
            return 0;
        }
        if (start <= end)
        {
            int middle = (start + end) / 2;
            if (uniqueIDs[middle].readNameID == readNameID)
            {
                return middle;
            }
            if (uniqueIDs[middle].readNameID > readNameID)
            {
                return searchReadNameID(readNameID, start, middle - 1);
            }
            return searchReadNameID(readNameID, middle + 1, end);
        }
        return start; // return the bigger one
    }

    /**
     * with a input readNameID, add it into readNameIDs.
     * if the input readNameID already exist in readNameIDs, return false. //这里修改了convert和unconvert
     */
    bool appendReadNameID(PosQuality &InBase, Alignment &InAlignment)
    {
        int idCount = uniqueIDs.size();
        if (idCount == 0 || InAlignment.readNameID > uniqueIDs.back().readNameID)
        {
            uniqueIDs.emplace_back(InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
        int index = searchReadNameID(InAlignment.readNameID, 0, idCount);
        if (uniqueIDs[index].readNameID == InAlignment.readNameID)
        {
            // if the new base is consistent with exist base's conversion status, ignore
            // otherwise, delete the exist conversion status
            if (uniqueIDs[index].removed)
            {
                return false;
            }
            if (uniqueIDs[index].isConverted != InBase.converted)
            {
                uniqueIDs[index].removed = true;
                if (uniqueIDs[index].isConverted)
                {
                    for (int i = 0; i < convertedQualities.size(); i++)
                    {
                        if (convertedQualities[i] == InBase.qual)
                        {
                            convertedQualities.erase(convertedQualities.begin() + i);
                            return false;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < unconvertedQualities.size(); i++)
                    {
                        if (unconvertedQualities[i] == InBase.qual)
                        {
                            unconvertedQualities.erase(unconvertedQualities.begin() + i);
                            return false;
                        }
                    }
                }
            }
            return false;
        }
        else
        {
            uniqueIDs.emplace(uniqueIDs.begin() + index, InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
    }

    /**
     * append the SAM information into this position. //这里修改了convert和unconvert
     */
    void appendBase(PosQuality &input, Alignment &a)
    {
        if (appendReadNameID(input, a))
        {
            if (input.converted)
            {
                convertedQualities += input.qual;
            }
            else
            {
                unconvertedQualities += input.qual;
            }
        }
    }
};
class Positions
{
public:
    ifstream refFile;
    // 全局变量
    ChromosomeFilePositions chromosomePos; // 存储染色体名称及其在文件中的流位置。用于快速在文件中查找新的染色体。
    long long int refCoveredPosition;      // 这是我们加载到refPositions中的参考染色体的最后位置。
    long long int location;                // 当前在参考染色体中的位置
    ostream *out_ = &cout;
    vector<Position *> refPositions;
    long long int reloadPos;
    string chromosome = " "; //正在使用的染色体块
    long long int loadingBlockSize = 1000000;


    /**
     * give a SAM line, extract the chromosome and position information.
     * return true if the SAM line is mapped. return false if SAM line is not maped.
     */
    bool getSAMChromosomePos(string *line, string &chr, long long int &pos)
    {
        int startPosition = 0;
        int endPosition = 0;
        int count = 0;

        while ((endPosition = line->find("\t", startPosition)) != string::npos)
        {
            if (count == 2)
            {
                chr = line->substr(startPosition, endPosition - startPosition);
            }
            else if (count == 3)
            {
                pos = stoll(line->substr(startPosition, endPosition - startPosition));
                if (chr == "*")
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }
            startPosition = endPosition + 1;
            count++;
        }
        return false;
    }

    bool getBamChromosomePos(bam1_t *b, bam_hdr_t *header, string &chr, long long int &pos)
    {
        // 检查是否未比对（FLAG 包含 BAM_FUNMAP）
        if (b->core.flag & BAM_FUNMAP)
        {
            chr = "*";
            pos = 0;
            return false;
        }

        // 获取染色体名称（需要 header 的 tid2name 映射）
        chr = string(header->target_name[b->core.tid]); // tid 是染色体ID

        // 获取比对位置（注意 BAM 是 1-based 坐标）
        pos = b->core.pos + 1; // bam1_core_t::pos 是 0-based，转为 1-based

        return true;
    }

    /**
     * given reference line (start with '>'), extract the chromosome information.
     * this is important when there is space in chromosome name. the SAM information only contain the first word.
     */
    string getChrName(string &inputLine)
    {
        string name;
        for (int i = 1; i < inputLine.size(); i++)
        {
            char c = inputLine[i];
            if (isspace(c))
            {
                break;
            }
            name += c;
        }

        // if (removedChrName)
        // {
        //     if (name.find("chr") == 0)
        //     {
        //         name = name.substr(3);
        //     }
        // }
        // else if (addedChrName)
        // {
        //     if (name.find("chr") != 0)
        //     {
        //         name = string("chr") + name;
        //     }
        // }
        return name;
    }

    /**
     * Scan the reference file. Record each chromosome and its position in file.
     */
    void LoadChromosomeNamesPos()
    {
        string line;
        while (refFile.good())
        {
            getline(refFile, line);
            if (line.front() == '>')
            { // this line is chromosome name
                chromosome = getChrName(line);
                streampos currentPos = refFile.tellg();
                chromosomePos.append(chromosome, currentPos);
            }
        }
        chromosomePos.sort();
        chromosome.clear();
    }

    void outputFunction(Position *pos)
    {
        *out_ << pos->Pchromosome << '\t'
              << to_string(pos->Plocation) << '\t'
              << pos->strand << '\t'
              << to_string(pos->convertedQualities.size()) << '\t'
              << to_string(pos->unconvertedQualities.size()) << '\n';
    }
    long long int samPos; // the position of current SAM line.
    /**
     * move the position which position smaller than refCoveredPosition - loadingBlockSize, output it.
     */
    void moveBlockToOutput()
    {
        if (refPositions.empty())
        {
            return;
        }
        int index;
        for (index = 0; index < refPositions.size(); index++)
        {
            if (refPositions[index]->Plocation < refCoveredPosition - loadingBlockSize)
            {
                if (refPositions[index]->empty() || refPositions[index]->strand == '?')
                {
                    // returnPosition(refPositions[index]);
                    delete refPositions[index];
                }
                else
                {
                    outputFunction(refPositions[index]);
                    delete refPositions[index]; // 节省内存，多线程不需要
                }
            }
            else
            {
                break;
            }
        }
        if (index != 0)
        {
            refPositions.erase(refPositions.begin(), refPositions.begin() + index);
        }
        int cmp = refPositions.size();
    }

    /**
     * move all the refPosition into output pool.
     */
    void moveAllToOutput()
    {
        if (refPositions.empty())
        {
            return;
        }
        for (int index = 0; index < refPositions.size(); index++)
        {
            if (refPositions[index]->empty() || refPositions[index]->strand == '?')
            {
                // returnPosition(refPositions[index]);
                delete refPositions[index];
            }
            else
            {
                vector<uniqueID>().swap(refPositions[index]->uniqueIDs);
                outputFunction(refPositions[index]);
                delete refPositions[index];
            }
        }
        refPositions.clear();
    }

    int getIndex(long long int &targetPos)
    {
        int firstPos = refPositions[0]->Plocation;
        return targetPos - firstPos;
    }

    /**
     * get a fasta line (not header), append the bases to positions.
     */
    void appendRefPosition(string &line)
    {
        char *b;
        for (int i = 0; i < line.size(); i++)
        {
            Position *newPos = new Position;
            newPos->set(chromosome, location + i);
            b = &line[i];
            if (*b == convertFrom)
            {
                newPos->set('+');
            }
            else if (*b == convertFromComplement)
            {
                newPos->set('-');
            }
            refPositions.push_back(newPos);
        }
        location += line.size();
    }

    void loadNewChromosome(string targetChromosome)
    {
        refFile.clear();
        // find the start position in file based on chromosome name.
        streampos startPos = chromosomePos.getChromosomePosInRefFile(targetChromosome);
        chromosome = targetChromosome;
        refFile.seekg(startPos, ios::beg);
        refCoveredPosition = 2 * loadingBlockSize;
        string line;
        location = 0;
        while (refFile.good())
        {
            getline(refFile, line);
            if (line.front() == '>')
            {           // this line is chromosome name
                return; // meet next chromosome, return it.
            }
            else
            {
                if (line.empty())
                {
                    continue;
                }
                appendRefPosition(line);
                if (location >= refCoveredPosition)
                {
                    return;
                }
            }
        }
    }

    void loadMore()
    {
        refCoveredPosition += loadingBlockSize;
        string line;
        while (refFile.good())
        {
            getline(refFile, line);
            if (line.front() == '>')
            { // meet next chromosome, return.
                return;
            }
            else
            {
                if (line.empty())
                {
                    continue;
                }
                appendRefPosition(line);
                if (location >= refCoveredPosition)
                {
                    return;
                }
            }
        }
    }

    /**
     * add position information from Alignment into ref position.
     */
    void appendPositions(Alignment &newAlignment)
    {
        if (!newAlignment.mapped || newAlignment.bases.empty())
        {
            return;
        }
        long long int startPos = newAlignment.location; // 1-based position
        // find the first reference position in pool.
        int index = getIndex(newAlignment.location);

        for (int i = 0; i < newAlignment.sequence.size(); i++)
        {
            PosQuality *b = &newAlignment.bases[i];
            if (b->remove)
            {
                continue;
            }
            Position *pos = refPositions[index + b->refPos];
            // assert (pos->Plocation == startPos + b->refPos);

            pos->appendBase(newAlignment.bases[i], newAlignment);
        }
    }
};
