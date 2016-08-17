CXX = /usr/bin/g++
CXXFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3 -I..

#OPENBLAS_DIR=../../../hpc-libs/openblas-git-cortexa57
OPENBLAS_DIR=../../../hpc-libs/openblas-git-armv8
OPENBLAS_IFLAGS = -I$(OPENBLAS_DIR)/include
OPENBLAS_LFLAGS = $(OPENBLAS_DIR)/lib/libopenblas.a -fopenmp

BLIS_DIR=../../../hpc-libs/blis-git-armv8a
BLIS_IFLAGS= -I$(BLIS_DIR)/include
BLIS_LFLAGS = $(BLIS_DIR)/lib/libblis.a -fopenmp

ATLAS_DIR=../../../hpc-libs/atlas
ATLAS_IFLAGS= -I$(ATLAS_DIR)/include
ATLAS_LFLAGS = $(ATLAS_DIR)/lib/libcblas.a $(ATLAS_DIR)/lib/libatlas.a

PERF_TARGETS = test-blocksizes perf-blis perf-openblas perf-atlas test-gemm-bliskernel
BLOCKSIZE_TARGETS = test-blocksizes

test-gemm-bliskernel:
	$(CXX) $(CXXFLAGS) $(BLOCKDEFS) $(BLIS_IFLAGS) -fpermissive -o $@ perf-gemm-bliskernel.cpp

test-blocksizes:
	@$(CXX) $(CXXFLAGS) $(BLOCKDEFS) -o $@ test-blocksizes.cpp

perf-libs:
	$(MAKE) perf-blis
	$(MAKE) perf-openblas
	$(MAKE) perf-atlas

perf-blis:
	$(MAKE) perf-blis-s
	@./$@ | tee $@-s.$${HOSTNAME}.txt
	$(MAKE) perf-blis-d
	@./$@ | tee $@-d.$${HOSTNAME}.txt

perf-blis-s: perf-blis.cpp
	g++ $(CXXFLAGS) $(BLIS_IFLAGS) -DPRECISION=1 -o perf-blis perf-blis.cpp $(BLIS_LFLAGS)

perf-blis-d: perf-blis.cpp
	g++ $(CXXFLAGS) $(BLIS_IFLAGS) -DPRECISION=2 -o perf-blis perf-blis.cpp $(BLIS_LFLAGS)

perf-atlas:
	$(MAKE) perf-atlas-s
	@./$@ | tee $@-s.$${HOSTNAME}.txt
	$(MAKE) perf-atlas-d
	@./$@ | tee $@-d.$${HOSTNAME}.txt

perf-atlas-s: perf-atlas.cpp
	g++ $(CXXFLAGS) $(ATLAS_IFLAGS) -DPRECISION=1 -o perf-atlas perf-atlas.cpp $(ATLAS_LFLAGS)

perf-atlas-d: perf-atlas.cpp
	g++ $(CXXFLAGS) $(ATLAS_IFLAGS) -DPRECISION=2 -o perf-atlas perf-atlas.cpp $(ATLAS_LFLAGS)

perf-openblas:
	$(MAKE) -e OPENBLAS_DIR=../../../hpc-libs/openblas-git-cortexa57 perf-openblas-s
	@./$@ | tee $@-a57-s.$${HOSTNAME}.txt
	$(MAKE) -e OPENBLAS_DIR=../../../hpc-libs/openblas-git-cortexa57 perf-openblas-d
	@./$@ | tee $@-a57-d.$${HOSTNAME}.txt
	$(MAKE) -e OPENBLAS_DIR=../../../hpc-libs/openblas-git-armv8 perf-openblas-s
	@./$@ | tee $@-generic-s.$${HOSTNAME}.txt
	$(MAKE) -e OPENBLAS_DIR=../../../hpc-libs/openblas-git-armv8 perf-openblas-d
	@./$@ | tee $@-generic-d.$${HOSTNAME}.txt

perf-openblas-s: perf-openblas.cpp
	g++ $(CXXFLAGS) $(OPENBLAS_IFLAGS) -DPRECISION=1 -o perf-openblas \
		perf-openblas.cpp $(OPENBLAS_LFLAGS)

perf-openblas-d: perf-openblas.cpp
	g++ $(CXXFLAGS) $(OPENBLAS_IFLAGS) -DPRECISION=2 -o perf-openblas \
	perf-openblas.cpp $(OPENBLAS_LFLAGS)

.PHONY: clean clean-perfdata runperf header $(PERF_TARGETS) $(BLOCKSIZE_TARGETS)

all:

header:
	@for i in $(BLOCKSIZE_TARGETS); \
	do echo "Prec   MC    KC    NC    MR    NR          blocked GEMM: t    MFLOPS" \
		| tee $$i.$${HOSTNAME}.$${PRECISION}.txt; \
	done

runperf: $(BLOCKSIZE_TARGETS)
	@for i in $(BLOCKSIZE_TARGETS); do ./$$i | tee -a $$i.$${HOSTNAME}.$${PRECISION}.txt; done

clean:
	rm -f $(PERF_TARGETS)

clean-perfdata:
	rm -f $(BLOCKSIZE_TARGETS:%=%.*.txt)

