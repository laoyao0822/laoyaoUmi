/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ALIGNMENT_3N_TABLE_H
#define ALIGNMENT_3N_TABLE_H

#include <string>
#include <sstream>
#include "utility_3n_table.h"
#include "sam.h"

extern char convertFrom;
extern char convertTo;
extern char convertFromComplement;
extern char convertToComplement;

using namespace std;

/**
 * the class to store information from one SAM line
 */
class Alignment
{
public:
    //string readName;
    string chromosome;
    long long int location;
    long long int mateLocation;
    int flag;
    bool mapped;
    char strand;
    string sequence;
    string quality;
    bool unique;
    string mapQ;
    int NH;
    int XM;
    int Zf;
    int Yf;
    vector<PosQuality> bases;
    CIGAR cigarString;
    MD_tag MD;
    unsigned long long readNameID;
    int sequenceCoveredLength; // the sum of number is cigarString;
    bool overlap;              // if the segment could overlap with the mate segment.
    bool paired;
    bool uniqueOnly;
    bool multipleOnly;

    void initialize()
    {
        chromosome.clear();
        location = -1;
        mateLocation = -1;
        flag = -1;
        mapped = false;
        MD.initialize();
        cigarString.initialize();
        sequence.clear();
        quality.clear();
        unique = false;
        mapQ.clear();
        NH = -1;
        bases.clear();
        readNameID = 0;
        sequenceCoveredLength = 0;
        overlap = false;
        paired = false;
    }

    /**
     * for start position in input Line, check if it contain the target information.
     */
    bool startWith(string *inputLine, int startPosition, string tag)
    {
        for (int i = 0; i < tag.size(); i++)
        {
            if (inputLine->at(startPosition + i) != tag[i])
            {
                return false;
            }
        }
        return true;
    }

    /**
     * generate a hash value for readName
     */
    void getNameHash(string &readName)
    {
        readNameID = 0;
        int a = 63689;
        for (size_t i = 0; i < readName.size(); i++)
        {
            readNameID = (readNameID * a) + (int)readName[i];
        }
    }

    /**
     * extract the information from SAM line to Alignment.
     */
    void parseInfo(string *line)
    {
        int startPosition = 0;
        int endPosition = 0;
        int count = 0;

        while ((endPosition = line->find("\t", startPosition)) != string::npos)
        {
            if (count == 0)
            {
                string readName = line->substr(startPosition, endPosition - startPosition);
                getNameHash(readName);
            }
            else if (count == 1)
            {
                flag = stoi(line->substr(startPosition, endPosition - startPosition));
                mapped = (flag & 4) == 0;
                paired = (flag & 1) != 0;
            }
            else if (count == 2)
            {
                chromosome = line->substr(startPosition, endPosition - startPosition);
            }
            else if (count == 3)
            {
                location = stoll(line->substr(startPosition, endPosition - startPosition));
            }
            else if (count == 4)
            {
                mapQ = line->substr(startPosition, endPosition - startPosition);
                if (mapQ == "1")
                {
                    unique = false;
                }
                else
                {
                    unique = true;
                }
            }
            else if (count == 5)
            {
                cigarString.loadString(line->substr(startPosition, endPosition - startPosition));
            }
            else if (count == 7)
            {
                mateLocation = stoll(line->substr(startPosition, endPosition - startPosition));
            }
            else if (count == 9)
            {
                sequence = line->substr(startPosition, endPosition - startPosition);
            }
            else if (count == 10)
            {
                quality = line->substr(startPosition, endPosition - startPosition);
            }
            else if (count > 10)
            {
                // YZ:A:-Yf:i:13  Zf:i:0  ZS:i:0  XN:i:0  XO:i:0  XG:i:0
                if (startWith(line, startPosition, "MD"))
                {
                    MD.loadString(line->substr(startPosition + 5, endPosition - startPosition - 5));
                }
                else if (startWith(line, startPosition, "NM"))
                {
                    NH = stoi(line->substr(startPosition + 5, endPosition - startPosition - 5));
                }
                else if (startWith(line, startPosition, "YZ"))
                {
                    strand = line->at(endPosition - 1);
                }
                else if (startWith(line, startPosition, "XM"))
                {
                    XM = stoi(line->substr(startPosition + 5, endPosition - startPosition - 5));
                }
                else if (startWith(line, startPosition, "Zf"))
                {
                    Zf = stoi(line->substr(startPosition + 5, endPosition - startPosition - 5));
                }
                else if (startWith(line, startPosition, "Yf"))
                {
                    Yf = stoi(line->substr(startPosition + 5, endPosition - startPosition - 5));
                }
            }
            startPosition = endPosition + 1;
            count++;
        }
        if (startWith(line, startPosition, "MD"))
        {
            MD.loadString(line->substr(startPosition + 5, endPosition - startPosition - 5));
        }
        else if (startWith(line, startPosition, "NM"))
        {
            NH = stoi(line->substr(startPosition + 5, endPosition - startPosition - 5));
        }
        else if (startWith(line, startPosition, "YZ"))
        {
            strand = line->at(endPosition - 1);
        }
    }

    void parseBamFields(bam1_t *b)
    {
        // 1. 提取 MD 字段（字符串）
        uint8_t *md_ptr = bam_aux_get(b, "MD");
        if (md_ptr)
        {
            MD.loadString(bam_aux2Z(md_ptr)); // MD:Z:6G12G0 -> "6G12G0"
        }
        // 2. 提取 NM 字段（整数）
        uint8_t *nm_ptr = bam_aux_get(b, "NM");
        if (nm_ptr)
        {
            NH = bam_aux2i(nm_ptr); // NM:i:2 -> 2
        }

        // 方法2：从 YZ:A:- 字段提取（如果存在）
        uint8_t *yz_ptr = bam_aux_get(b, "YZ");
        if (yz_ptr && yz_ptr[0] == 'A')
        {
            strand = bam_aux2A(yz_ptr); // YZ:A:- -> '-'
        }

        // 4. 提取 XM（整数）
        uint8_t *xm_ptr = bam_aux_get(b, "XM");
        if (xm_ptr)
        {
            XM = bam_aux2i(xm_ptr); // XM:i:1 -> 1
        }

        // 5. 提取 Zf（整数）
        uint8_t *zf_ptr = bam_aux_get(b, "Zf");
        if (zf_ptr)
        {
            Zf = bam_aux2i(zf_ptr); // Zf:i:0 -> 0
        }

        // 6. 提取 Yf（整数）
        uint8_t *yf_ptr = bam_aux_get(b, "Yf");
        if (yf_ptr)
        {
            Yf = bam_aux2i(yf_ptr); // Yf:i:13 -> 13
        }
    }

    void parseInfo(bam1_t *b, bam_hdr_t *header)
    {
        // 初始化所有成员变量
        initialize();

        // 检查输入有效性
        if (b == nullptr || b->data == nullptr)
        {
            throw std::runtime_error("Invalid bam1_t structure");
        }

        // 1. 获取查询名称 (QNAME)
        string readName = bam_get_qname(b);
        if (readName.empty())
        {
            throw std::runtime_error("Empty read name");
        }
        getNameHash(readName);

        // 2. 解析 FLAG
        flag = b->core.flag;
        mapped = (flag & BAM_FUNMAP) == 0;
        paired = (flag & BAM_FPAIRED) != 0;

        // 3. 获取参考序列名称 (RNAME)
        if (b->core.tid >= 0 && header && header->target_name)
        {
            chromosome = header->target_name[b->core.tid];
        }
        else
        {
            chromosome = "*";
        }

        // 4. 获取位置 (POS) - 转换为1-based
        location = b->core.pos + 1;

        // 5. 获取 mapping quality (MAPQ)
        mapQ = std::to_string(b->core.qual);
        if (mapQ == "1")
        {
            unique = false;
        }
        else
        {
            unique = true;
        }

        // 6. 解析 CIGAR 字符串
        uint32_t *cigar = bam_get_cigar(b);
        if (cigar && b->core.n_cigar > 0)
        {
            std::stringstream cigarSS;
            for (int i = 0; i < b->core.n_cigar; ++i)
            {
                int op = bam_cigar_op(cigar[i]);
                int len = bam_cigar_oplen(cigar[i]);
                cigarSS << len << "MIDNSHP=XB"[op];
            }
            cigarString.loadString(cigarSS.str());
        }
        else
        {
            cigarString.loadString("*");
        }

        // 7. 获取 mate 位置 (PNEXT)
        mateLocation = b->core.mpos + 1;

        // 8. 获取序列 (SEQ)
        uint8_t *seq = bam_get_seq(b);
        sequence.clear();
        if (seq && b->core.l_qseq > 0)
        {
            for (int i = 0; i < b->core.l_qseq; ++i)
            {
                sequence += "=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)];
            }
        }

        // 9. 获取质量值 (QUAL)
        uint8_t *qual = bam_get_qual(b);
        quality.clear();
        if (qual && b->core.l_qseq > 0)
        {
            for (int i = 0; i < b->core.l_qseq; ++i)
            {
                quality += (qual[i] + 33); // 转换为ASCII
            }
        }

        // 10.可选字段
        parseBamFields(b);

        // 11. 获取链方向 (从flag判断)
        strand = (flag & BAM_FREVERSE) ? '-' : '+';

        // 根据uniqueOnly/multipleOnly条件过滤
        if ((uniqueOnly && !unique) || (multipleOnly && unique))
        {
            return;
        }

        appendBase();
    }

    /**
     * change the overlap = true, if the read is not uniquely mapped or the read segment is overlap to it's mate.
     */
    void checkOverlap()
    {
        if (!unique)
        {
            overlap = true;
        }
        else
        {
            if (paired && (location + sequenceCoveredLength >= mateLocation))
            {
                overlap = true;
            }
            else
            {
                overlap = false;
            }
        }
    }

    /**
     * parse the sam line to alignment information
     */
    void parse(string *line)
    {
        initialize();
        parseInfo(line);
        if ((uniqueOnly && !unique) || (multipleOnly && unique))
        {
            return;
        }
        appendBase();
    }

    /**
     *  scan all base in read sequence label them if they are qualified.
     */
    void appendBase()
    {
        if (!mapped || sequenceCoveredLength > 500000)
        { // if the read's intron longer than 500,000 ignore this read
            return;
        }

        bases.reserve(sequence.size());
        for (int i = 0; i < sequence.size(); i++)
        {
            bases.emplace_back(i);
        }
        int pos = adjustPos();

        string match;
        while (MD.getNextSegment(match))
        {
            if (isdigit(match.front()))
            { // the first char of match is digit this is match
                int len = stoi(match);
                for (int i = 0; i < len; i++)
                {
                    while (bases[pos].remove)
                    {
                        pos++;
                    }
                    if ((strand == '+' && sequence[pos] == convertFrom) ||
                        (strand == '-' && sequence[pos] == convertFromComplement))
                    {
                        bases[pos].setQual(quality[pos], false);
                    }
                    else
                    {
                        bases[pos].remove = true;
                    }
                    pos++;
                }
            }
            else if (isalpha(match.front()))
            { // this is mismatch or conversion
                char refBase = match.front();
                // for + strand, it should have C->T change
                // for - strand, it should have G->A change
                while (bases[pos].remove)
                {
                    pos++;
                }

                if ((strand == '+' && refBase == convertFrom && sequence[pos] == convertTo) ||
                    (strand == '-' && refBase == convertFromComplement && sequence[pos] == convertToComplement))
                {
                    bases[pos].setQual(quality[pos], true);
                }
                else
                {
                    bases[pos].remove = true;
                }
                pos++;
            }
            else
            { // deletion. do nothing.
            }
        }
    }

    /**
     * adjust the reference position in bases
     */
    int adjustPos()
    {

        int readPos = 0;
        int returnPos = 0;
        int seqLength = sequence.size();

        char cigarSymbol;
        int cigarLen;
        sequenceCoveredLength = 0;

        while (cigarString.getNextSegment(cigarLen, cigarSymbol))
        {
            sequenceCoveredLength += cigarLen;
            if (cigarSymbol == 'S')
            {
                if (readPos == 0)
                { // soft clip is at the begin of the read
                    returnPos = cigarLen;
                    for (int i = cigarLen; i < seqLength; i++)
                    {
                        bases[i].refPos -= cigarLen;
                    }
                }
                else
                { // soft clip is at the end of the read
                  // do nothing
                }
                readPos += cigarLen;
            }
            else if (cigarSymbol == 'N')
            {
                for (int i = readPos; i < seqLength; i++)
                {
                    bases[i].refPos += cigarLen;
                }
            }
            else if (cigarSymbol == 'M')
            {
                for (int i = readPos; i < readPos + cigarLen; i++)
                {
                    bases[i].remove = false;
                }
                readPos += cigarLen;
            }
            else if (cigarSymbol == 'I')
            {
                for (int i = readPos + cigarLen; i < seqLength; i++)
                {
                    bases[i].refPos -= cigarLen;
                }
                readPos += cigarLen;
            }
            else if (cigarSymbol == 'D')
            {
                for (int i = readPos; i < seqLength; i++)
                {
                    bases[i].refPos += cigarLen;
                }
            }
        }
        return returnPos;
    }

    bool parse_cigar(int &qlen, int &sclen)
    {
        qlen = 0;
        sclen = 0;

        int num = 0;
        for (char c : cigarString.s)
        {
            if (isdigit(c))
            {
                num = num * 10 + (c - '0');
            }
            else
            {
                if (c == 'M' || c == 'I' || c == 'S')
                    qlen += num;
                if (c == 'S')
                    sclen += num;
                num = 0;
            }
        }
        return (XM * 20 <= (qlen - sclen)) && (Zf <= 3) && (3 * Zf <= Zf + Yf);
    }
};

#endif // ALIGNMENT_3N_TABLE_H
