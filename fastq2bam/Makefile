#include ../makefile.defs


CXX      = g++ #-std=c++11 
BAMTOOLS = $(realpath ../bamtools/)
LIBGAB   = $(realpath ../libgab/)

CXXFLAGS  = -lm -O3 -Wall  -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${LIBGAB}VCFparser/  -I${BAMTOOLS}/src/ -c
LDFLAGS  = ${BAMTOOLS}/build/src/api/libbamtools.a -lz


all: fastq2bam

fastq2bam.o:	fastq2bam.cpp
	${CXX} ${CXXFLAGS} fastq2bam.cpp

fastq2bam:	fastq2bam.o ${LIBGAB}/libgab.a ${LIBGAB}/PutProgramInHeader.o ${LIBGAB}/gzstream/gzstream.o 
	${CXX} -o $@ $^ $(LDFLAGS)


clean :
	rm -f *.o fastq2bam

