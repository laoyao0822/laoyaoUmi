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

#include <iostream>
#include <getopt.h>
// #include "position_3n_table.h"
#include "alignment_3n_table.h"
#include "position_3n_table.h"
#include <vector>
#include <fstream>
#include <cassert>
#include <cstdio>
#include <cstring>
#include "sam.h"
using namespace std;

string alignmentFileName;
bool standardInMode = false;
string refFileName;
char convertFrom = 'C';
char convertTo = 'T';
char convertFromComplement = asc2dnacomp[convertFrom];
char convertToComplement = asc2dnacomp[convertTo];
long long int loadingBlockSize = 1000000;
extern const unsigned int CHUNK_SIZE;

extern sam_hdr_t *bam_header;
extern bam1_t ***b_array;

extern bam1_t *getBam1_t(unsigned int index);




void hisat_3n_table(unsigned int *merge, int laoyaoStart, int laoyaoEnd,
                   string outputFileName,
                   bool uniqueOnly, bool multipleOnly, bool needFilter)
{
    Positions positions;
    ofstream tableFile;
    if (!outputFileName.empty())
    {
        tableFile.open(outputFileName, ios_base::out);
        positions.out_ = &tableFile;
    } else {
        cout << "please check the outputFile Path for " << outputFileName << endl; 
        return;
    }

    positions.refFile.open("/root/m5C/data/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa", ios_base::in); //硬编码！请注意参考文件路径
    positions.LoadChromosomeNamesPos();

    string chromosome;
    int qlen, sclen;
    long long int samPos;


    bam1_t *b = bam_init1();

    Alignment *newAlignment = new Alignment();
    newAlignment->uniqueOnly = uniqueOnly;
    newAlignment->multipleOnly = multipleOnly;


    for (int i = laoyaoStart; i <= laoyaoEnd; i++)
    {
        b = getBam1_t(merge[i]);
        if (!positions.getBamChromosomePos(b, bam_header, chromosome, samPos))
        {
            continue;
        }
        if (positions.chromosome != chromosome)
        {
            positions.moveAllToOutput();
            positions.loadNewChromosome(chromosome);
            positions.reloadPos = loadingBlockSize;
        }
        while (samPos > positions.reloadPos)
        {
            positions.moveBlockToOutput();
            positions.loadMore();
            positions.reloadPos += loadingBlockSize;
        }
        newAlignment->parseInfo(b, bam_header);
        if (needFilter && !newAlignment->parse_cigar(qlen, sclen))
        {
            continue;
        }
        positions.appendPositions(*newAlignment);
    }

    tableFile.close();
}
