#include ../makefile.defs


CXX      = g++ #-std=c++11 
BAMTOOLS = $(realpath ../bamtools/)
LIBGAB   = $(realpath ../libgab/)

CXXFLAGS  = -lm -O3 -Wall  -I${LIBGAB} -I${LIBGAB}VCFparser/ -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/ -c
LDFLAGS  = ${BAMTOOLS}/lib/libbamtools.a -lz


all: fastq2bam

fastq2bam.o:	fastq2bam.cpp
	${CXX} ${CXXFLAGS} fastq2bam.cpp

fastq2bam:	fastq2bam.o ${LIBGAB}/utils.o ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/gzstream/gzstream.o ${LIBGAB}/FastQObj.o  ${LIBGAB}/FastQParser.o
	${CXX} -o $@ $^ $(LDFLAGS)


clean :
	rm -f *.o fastq2bam

