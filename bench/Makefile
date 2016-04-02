CPP = /usr/bin/g++
CPPFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3 -I..

OPENBLAS_DIR=../../../hpc-libs/openblas
OPENBLAS_IFLAGS = -I$(OPENBLAS_DIR)/include
OPENBLAS_LFLAGS = $(OPENBLAS_DIR)/lib/libopenblas.a -fopenmp

PERF_TARGETS = test-blocksizes

perf-blis: perf-blis.cpp
	g++  -Wall -Wextra -Wpedantic -std=c++11  -I/home/chris/workspace/bachelorarbeit/hpc-libs/blis/include/blis/ perf-blis.cpp   -o perf-blis temp/libblis.a -fopenmp

test-blocksizes:
	g++ $(CPPFLAGS) $(BLOCKDEFS) -o test-blocksizes test-blocksizes.cpp

perf-openblas: perf-openblas.cpp
	g++ $(CPPFLAGS) $(OPENBLAS_IFLAGS) -o perf-openblas perf-openblas.cpp $(OPENBLAS_LFLAGS)

.PHONY: clean clean-perfdata runperf header $(PERF_TARGETS)

ALL:

header:
	for i in $(PERF_TARGETS); do echo "   MC    KC    NC    MR    NR     blocked GEMM: t    MFLOPS" | tee $$i.$${HOSTNAME}.out; done

runperf: $(PERF_TARGETS)
	for i in $(PERF_TARGETS); do ./$$i | tee -a $$i.$${HOSTNAME}.out; done

clean:
	rm -f $(PERF_TARGETS)

clean-perfdata:
	rm -f $(PERF_TARGETS:%=%.*.out)

