#include ../makefile.defs

CXX      = g++
BAMTOOLS = $(realpath ../bamtools/)
LIBGAB   = $(realpath ../libgab/)



#LIBGAB   = /home/gabriel_renaud/lib/
LDLIBS+= ${BAMTOOLS}/build/src/api/libbamtools.a ${LIBGAB}/libgab.a -lm -lz

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/src/  -c
#LDFLAGS  = -lz


all: bam2fastq 

bam2fastq.o:	bam2fastq.cpp
	${CXX} ${CXXFLAGS} bam2fastq.cpp


bam2fastq:	bam2fastq.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f bam2fastq.o bam2fastq

