# for static compilation: STATIC=1 make
# for debug: make DEBUG=1
# dynamic linking requires module load gcc/4.9.2 samtools gsl atlas lapack armadillo
# static linking requires module load gcc/4.9.2 samtools gsl atlas lapack armadillo zlib boost

SUFFIX = # set on command line for static compilation
MACROS = -DVERSION=$(shell echo `cat VERSION`)
ifdef DEBUG
  MACROS += -DDEBUG
endif

BINDIR=../bin

CC = gcc
CXX = g++
LIBS = -lgsl -lgslcblas -lboost_regex -lhts -lboost_iostreams -lz -llapack -lblas -lgfortran  
OMP = -fopenmp
OPT = -O3
OPTSAFE = -O2
#CXXFLAGS = -std=c++11

ifeq ($(shell echo `uname -s`),Darwin)
#  CLEANSUFFIX = $(shell echo `cat VERSION`-`uname -s`-`ar -x /opt/local/lib/libgsl.a rng.o && file rng.o | sed -E 's/.*\ (.*)/\1/'`)
   CLEANSUFFIX = -mac # $(shell echo -`cat VERSION`-`uname -s`-`uname -m`)
endif
ifeq ($(shell echo `uname -s`),Linux)
  CLEANSUFFIX = -linux # $(shell echo -`cat VERSION`-`uname -s`-`uname -m`)
endif
ifdef STATIC
  CC += -static-libgcc 
  CXX += -static-libgcc 
  SUFFIX = ${CLEANSUFFIX}
  ifeq ($(shell echo `uname -s`),Darwin)
    $(warning "Static compilation for Mac OS X")
# comment out #define ARMA_USE_WRAPPER in /usr/local/include/armadillo_bits/config.hpp
    CXXFLAGS += -framework Accelerate -mmacosx-version-min=10.6
    LIBS = /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a /usr/local/lib/libboost_iostreams.a -lz /usr/local/lib/libboost_regex-mt.a /usr/local/lib/libbam.a 
  endif
  ifeq ($(shell echo `uname -s`),Linux)
    $(warning "Static compilation for Linux")
    CXXFLAGS += -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
		LIBS += -lquadmath
  endif
endif
  

.PHONY : all
all : mmseq bam2hits extract_transcripts offsetGTF hitstools t2g_hits mmcollapse mmdiff

mmseq : sokal.o uh.o hitsio.o mmseq.cpp
	$(CXX) ${MACROS} $(CXXFLAGS) ${OMP} ${OPT} $^ -o ${BINDIR}/$@${SUFFIX} ${LIBS}

bam2hits : fasta.o hitsio.o bam2hits.cpp
	${CXX} ${MACROS} $(CXXFLAGS) ${OMP} ${OPT} $^ -o ${BINDIR}/$@${SUFFIX} ${LIBS}

extract_transcripts : fasta.o extract_transcripts.cpp
	${CXX} ${MACROS} $(CXXFLAGS) ${OPT} $^ -o ${BINDIR}/$@${SUFFIX}

offsetGTF : offsetGTF.cpp
	${CXX} ${MACROS} $(CXXFLAGS) ${OPT} $^ -o ${BINDIR}/$@${SUFFIX}

hitstools: hitsio.o hitstools.cpp
	$(CXX) ${MACROS} $(CXXFLAGS) ${OMP} ${OPT} $^ -o ${BINDIR}/$@${SUFFIX} ${LIBS}

t2g_hits : t2g_hits.cpp
	${CXX} ${MACROS} $(CXXFLAGS) ${OPT} $^ -o ${BINDIR}/$@${SUFFIX}

mmcollapse: uh.o sokal.o mmcollapse.cpp
	${CXX} ${MACROS} -DARMA_64BIT_WORD ${CXXFLAGS} ${OMP} ${OPT} $^ -o ${BINDIR}/$@${SUFFIX} ${LIBS}

mmdiff: bms.o mmdiff.cpp
	${CXX} ${MACROS} ${CXXFLAGS} ${OMP} ${OPTSAFE} $^ -o ${BINDIR}/$@${SUFFIX} ${LIBS}

sokal.o : sokal.cc sokal.hh
	$(CXX) ${OPT} -c sokal.cc

fasta.o : fasta.c fasta.h

uh.o : uh.cpp uh.hh
	$(CXX) ${OMP} ${OPT} -c uh.cpp 

hitsio.o: hitsio.cpp hitsio.hpp
	$(CXX) ${MACROS} $(CXXFLAGS) ${OMP} ${OPT} $^ -c hitsio.cpp

bms.o: bms.cpp bms.hpp
	$(CXX) ${MACROS} $(CXXFLAGS) ${OMP} ${OPTSAFE} $^ -c bms.cpp

.PHONY : clean

clean :
	rm -f ${BINDIR}/mmseq ${BINDIR}/bam2hits ${BINDIR}/extract_transcripts ${BINDIR}/hitstools ${BINDIR}/t2g_hits ${BINDIR}/mmcollapse ${BINDIR}/mmseq${CLEANSUFFIX} ${BINDIR}/bam2hits${CLEANSUFFIX} ${BINDIR}/extract_transcripts${CLEANSUFFIX} ${BINDIR}/hitstools${CLEANSUFFIX} ${BINDIR}/t2g_hits${CLEANSUFFIX} ${BINDIR}/mmcollapse${CLEANSUFFIX} ${BINDIR}/mmdiff${CLEANSUFFIX} *.o *.hpp.gch
