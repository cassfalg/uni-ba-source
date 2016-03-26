CPP = /usr/bin/g++
CPPFLAGS = -Wall -Wextra -Wpedantic -std=c++11

PERF_TARGETS = linear-chase fused-linear-chase random-chase

.PHONY: clean clean-perfdata runperf

ALL:

runperf: $(PERF_TARGETS)
	for i in $(PERF_TARGETS); do ./$$i | tee $$i.$${HOSTNAME}.out; done

clean:
	rm -f $(PERF_TARGETS)

clean-perfdata:
	rm -f $(PERF_TARGETS:%=%.*.out)