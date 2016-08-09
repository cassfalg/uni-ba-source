CXX = /usr/bin/g++
CXXFLAGS += -Wall -Wextra -Wpedantic -std=c++11 -O3 -I..

OPENBLAS_DIR=../../../hpc-libs/openblas-git-armv8
OPENBLAS_IFLAGS = -I$(OPENBLAS_DIR)/include
OPENBLAS_LFLAGS = $(OPENBLAS_DIR)/lib/libopenblas.a -fopenmp

BLIS_DIR=../../../hpc-libs/blis-git-armv8a
BLIS_IFLAGS= -I$(BLIS_DIR)/include
BLIS_LFLAGS = $(BLIS_DIR)/lib/libblis.a -fopenmp

ATLAS_DIR=../../../hpc-libs/atlas
ATLAS_IFLAGS= -I$(ATLAS_DIR)/include
ATLAS_LFLAGS = $(ATLAS_DIR)/lib/libcblas.a $(ATLAS_DIR)/lib/libatlas.a

all: tests

.PHONY: clean tests all compile

MY_TARGETS = refblas_vs_myblas_refugemm

TARGETS = $(MY_TARGETS) openblas_vs_myblas_refugemm \
			blis_vs_myblas_refmygemm atlas_vs_myblas_refugemm

atlas_vs_myblas_refugemm: atlas_vs_myblas_refugemm.cpp
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) $(ATLAS_IFLAGS) -o $@ $@.cpp $(ATLAS_LFLAGS)

openblas_vs_myblas_refugemm: openblas_vs_myblas_refugemm.cpp
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) $(OPENBLAS_IFLAGS) -o $@ $@.cpp $(OPENBLAS_LFLAGS)

blis_vs_myblas_refmygemm: blis_vs_myblas_refmygemm.cpp
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) $(BLIS_IFLAGS) -o $@ $@.cpp $(BLIS_LFLAGS)

$(MY_TARGETS): $(MY_TARGETS).cpp
	@$(CXX) $(CXXFLAGS) $(PP_FLAGS) -o $@ $@.cpp

compile: $(TARGETS)

tests: $(TARGETS)
	for i in $(TARGETS); do ./$$i | tee $$i.$${HOSTNAME}.$${PRECISION}.txt; done

clean:
	rm -f $(TARGETS) *.txt

