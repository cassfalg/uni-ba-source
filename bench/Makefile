CPP = /usr/bin/g++
CPPFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3 -I..

OPENBLAS_DIR=../../../hpc-libs/openblas
OPENBLAS_IFLAGS = -I$(OPENBLAS_DIR)/include
OPENBLAS_LFLAGS = $(OPENBLAS_DIR)/lib/libopenblas.a -fopenmp

BLIS_DIR=../../../hpc-libs/blis
BLIS_IFLAGS= -I$(BLIS_DIR)/include
BLIS_LFLAGS = $(BLIS_DIR)/lib/libblis.a -fopenmp

ATLAS_DIR=../../../hpc-libs/atlas
ATLAS_IFLAGS= -I$(ATLAS_DIR)/include
ATLAS_LFLAGS = $(ATLAS_DIR)/lib/libcblas.a $(ATLAS_DIR)/lib/libatlas.a

PERF_TARGETS = test-blocksizes perf-blis perf-openblas perf-atlas
BLOCKSIZE_TARGETS = test-blocksizes

perf-blis-old: perf-blisold.cpp
	g++  -Wall -Wextra -Wpedantic -std=c++11  -I/home/chris/workspace/bachelorarbeit/hpc-libs/blis/include/blis/ perf-blis.cpp   -o perf-blis temp/libblis.a -fopenmp

test-blocksizes:
	g++ $(CPPFLAGS) $(BLOCKDEFS) -o test-blocksizes test-blocksizes.cpp

perf-blis: perf-blis.cpp
	g++ $(CPPFLAGS) $(BLIS_IFLAGS) -o perf-blis perf-blis.cpp $(BLIS_LFLAGS)

perf-atlas: perf-atlas.cpp
	g++ $(CPPFLAGS) $(ATLAS_IFLAGS) -o perf-atlas perf-atlas.cpp $(ATLAS_LFLAGS)

perf-openblas: perf-openblas.cpp
	g++ $(CPPFLAGS) $(OPENBLAS_IFLAGS) -o perf-openblas perf-openblas.cpp $(OPENBLAS_LFLAGS)

.PHONY: clean clean-perfdata runperf header $(PERF_TARGETS) $(BLOCKSIZE_TARGETS)

all:

header:
	for i in $(BLOCKSIZE_TARGETS); do echo "   MC    KC    NC    MR    NR     blocked GEMM: t    MFLOPS" | tee $$i.$${HOSTNAME}.out; done

runperf: $(BLOCKSIZE_TARGETS)
	for i in $(BLOCKSIZE_TARGETS); do ./$$i | tee -a $$i.$${HOSTNAME}.out; done

clean:
	rm -f $(PERF_TARGETS)

clean-perfdata:
	rm -f $(BLOCKSIZE_TARGETS:%=%.*.out)

