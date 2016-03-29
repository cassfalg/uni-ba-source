CPP = /usr/bin/g++
CPPFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3 


PERF_TARGETS = test-blocksizes

test-blocksizes:
	g++ $(CPPFLAGS) $(BLOCKDEFS) -o test-blocksizes test-blocksizes.cpp

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

