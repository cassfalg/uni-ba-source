CPP = /usr/bin/g++
CPPFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3 -I..


PERF_TARGETS = test-gemm


.PHONY: clean test

ALL: test

test: $(PERF_TARGETS)
	for i in $(PERF_TARGETS); do ./$$i; done

clean:
	rm -f $(PERF_TARGETS)

