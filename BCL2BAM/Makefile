#include ../makefile.defs

CXX      = g++
BAMTOOLS = $(realpath ../bamtools/)
LIBGAB   = $(realpath ../libgab/)

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/   -I${BAMTOOLS}/build/src/ -I${BAMTOOLS}/src/  -c
LDLIBS   += ${BAMTOOLS}/build/src/libbamtools.a -lm -lz


all: bcl2bam 

bcl2bam.o:	bcl2bam.cpp
	${CXX} ${CXXFLAGS} bcl2bam.cpp

bcl2bam:	bcl2bam.o ${LIBGAB}/libgab.a   ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS)

clean :
	rm -f bcl2bam.o bcl2bam

