CXX = /usr/bin/g++
CXXFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3 -I..
PRECISION ?= 2

OPENBLAS_DIR=../../../hpc-libs/openblas-git-armv8
#OPENBLAS_DIR=../../../hpc-libs/openblas-git-cortexa57/
OPENBLAS_IFLAGS = -I$(OPENBLAS_DIR)/include
OPENBLAS_LFLAGS = $(OPENBLAS_DIR)/lib/libopenblas.a -fopenmp

#BLIS_DIR=../../../hpc-libs/blis-git-armv8a
BLIS_DIR=/home/chris/workspace/bachelorarbeit/hpc-libs/blis-0.1.8-armv8a
BLIS_IFLAGS= -I$(BLIS_DIR)/include
BLIS_LFLAGS = $(BLIS_DIR)/lib/libblis.a -fopenmp

ATLAS_DIR=../../../hpc-libs/atlas
ATLAS_IFLAGS= -I$(ATLAS_DIR)/include
ATLAS_LFLAGS = $(ATLAS_DIR)/lib/libcblas.a $(ATLAS_DIR)/lib/libatlas.a

MY_TARGETS = refblas_vs_myblas_refugemm

TARGETS = $(MY_TARGETS) \
			openblas_vs_myblas_refugemm \
			openblas_armv8_vs_myblas_refugemm \
			openblas_cortexa57_vs_myblas_refugemm \
			blis_vs_myblas_refugemm \
			atlas_vs_myblas_refugemm \
			refblas_vs_myblas_simple_4x4_cm_ugemm \
			refblas_vs_myblas_simple_4x4_rm_ugemm

KERNELS = mark2a mark2b

all: tests

.PHONY: $(TARGETS) clean tests all compile manual_tests compile_kernels

atlas_vs_myblas_refugemm:
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) $(ATLAS_IFLAGS) -o $@ $@.cpp $(ATLAS_LFLAGS)
#	@./$@ | tee $@.$${HOSTNAME}.$(PRECISION).txt

openblas_vs_myblas_refugemm:
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) $(OPENBLAS_IFLAGS) -o $@ $@.cpp $(OPENBLAS_LFLAGS)
#	./$@ | tee $@.$${HOSTNAME}.$(PRECISION).txt

openblas_armv8_vs_myblas_refugemm:
	@$(MAKE) -e OPENBLAS_DIR=../../../hpc-libs/openblas-git-armv8 openblas_vs_myblas_refugemm
	@./openblas_vs_myblas_refugemm \
			| tee -a openblas_armv8_vs_myblas_refugemm.$${HOSTNAME}.$(PRECISION).txt

openblas_cortexa57_vs_myblas_refugemm:
	@$(MAKE) -e OPENBLAS_DIR=../../../hpc-libs/openblas-git-cortexa57 openblas_vs_myblas_refugemm
	@./openblas_vs_myblas_refugemm \
			| tee -a openblas_cortexa57_vs_myblas_refugemm.$${HOSTNAME}.$(PRECISION).txt

blis_vs_myblas_refugemm:
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) $(BLIS_IFLAGS) -o $@ $@.cpp $(BLIS_LFLAGS)
#	@./$@ | tee $@.$${HOSTNAME}.$(PRECISION).txt

$(MY_TARGETS): $(MY_TARGETS).cpp
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) -o $@ $@.cpp
#	@./$@ | tee $@.$${HOSTNAME}.$(PRECISION).txt

refblas_vs_myblas_simple_4x4_cm_ugemm: refblas_vs_myblas_simple_4x4_cm_ugemm.cpp \
				../hpc/myblas/kernels/my_dgemm_simple_4x4_cm.s \
				../hpc/myblas/kernels/my_dgemm_simple_4x4_cm.h
	@as -o my_dgemm_simple_4x4_cm.o ../hpc/myblas/kernels/my_dgemm_simple_4x4_cm.s
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) my_dgemm_simple_4x4_cm.o -o $@ $@.cpp
#	@./$@ | tee $@.$${HOSTNAME}.$(PRECISION).txt

refblas_vs_myblas_simple_4x4_rm_ugemm: refblas_vs_myblas_simple_4x4_rm_ugemm.cpp \
				../hpc/myblas/kernels/my_dgemm_simple_4x4_rm.s \
				../hpc/myblas/kernels/my_dgemm_simple_4x4_rm.h
	@as -o my_dgemm_simple_4x4_rm.o ../hpc/myblas/kernels/my_dgemm_simple_4x4_rm.s
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) my_dgemm_simple_4x4_rm.o -o $@ $@.cpp
#	@./$@ | tee $@.$${HOSTNAME}.$(PRECISION).txt

#refblas_vs_myblas_mark2a: refblas_vs_myblas_mark2a.cpp \
#				../hpc/myblas/kernels/my_dgemm_mark2a.s \
#				../hpc/myblas/kernels/my_dgemm_mark2a.h
#	@as -o my_dgemm_mark2a.o ../hpc/myblas/kernels/my_dgemm_mark2a.s
#	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) my_dgemm_mark2a.o -o $@ $@.cpp
#	@./$@ | tee $@.$${HOSTNAME}.$(PRECISION).txt

%.o : %.s #cancel default rule

mark%.o: ../hpc/myblas/kernels/mark%.s
	as -o $@ $<


refblas_vs_myblas_mark%: mark%.o
	$(CXX) $(CXXFLAGS) $(PP_FLAGS) $< -o $@ refblas_vs_myblas_kernels.cpp

test_kernels: refblas_vs_myblas_$(KERNELS)
	for i in $(KERNELS); do ./refblas_vs_myblas_$$i \
		| tee -a refblas_vs_myblas_$$i.$${HOSTNAME}.$(PRECISION).txt; done

compile: $(TARGETS)

tests: $(TARGETS)
	for i in $(TARGETS); do ./$$i | tee -a $$i.$${HOSTNAME}.$(PRECISION).txt; done

manual_tests:
	@as -o my_dgemm_simple_4x4_cm.o ../hpc/myblas/kernels/my_dgemm_simple_4x4_cm.s
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) my_dgemm_simple_4x4_cm.o -o $@ $@.cpp
	./$@

clean:
	rm -f $(TARGETS) refblas_vs_myblas_$(KERNELS) *.o

allclean: clean
	rm -f *.txt


