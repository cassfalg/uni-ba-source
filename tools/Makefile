CPP = /usr/bin/g++
CPPFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3

PERF_TARGETS = random-chase linear-chase fused-linear-chase

.PHONY: clean clean-perfdata runperf


ALL: runperf

runperf: $(PERF_TARGETS)
	for i in $(PERF_TARGETS); do ./$$i | tee $$i.$${HOSTNAME}.txt; done

clean:
	rm -f $(PERF_TARGETS)

clean-perfdata:
	rm -f $(PERF_TARGETS:%=%.*.out)

