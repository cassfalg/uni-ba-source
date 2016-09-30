#include <cstdio>
#include <hpc/util/walltime.h>

#ifndef MIN_SIZE
#define MIN_SIZE 20000u
#endif
#ifndef GRANULARITY
#define GRANULARITY 10000u
#endif
#ifndef MAX_SIZE
#define MAX_SIZE 500000u // uses about 1,5 GiB of memory, about all we have
#endif

extern "C" {
	void fmla_test(long num_iter);
}

int
main() {
    hpc::util::WallTime<double> wallTime;
    long count = 300000000;
    printf("%20s %10s %10s %10s\n", "Instruction", "Count", "Time", "IPC");
    wallTime.tic();
    fmla_test(count);
    double t = wallTime.toc();
    // Instructions per cycle: 2 Ghz equals 2 000 000 000 cycles per second
    // function does 30 (+ sub + branch) instructions per iteration
    printf("%20s %10ld %10.2lf %10.5lf\n", "fmla Vt, Vt, Vt"
    		, count, t, count/2.0/1000000000.0*18/t);
}
