INC = ~/

models = flo_evol.o flo_ising_random_simp.o flo_ising_random_simp_shift_real.o

small_methods = matrix_algebra.o screen_output.o mersenne.o

transitions = basis_transition.o basic_full.o basic_parity.o parity_full.o

test_objects = test.o $(models) $(small_methods) $(transitions)

CXXFLAGS = -O3 -fopenmp
CXX = g++

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
CXX = icpc
CXXFLAGS += $(INC)
endif

all: mbl

mbl: $(main_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

test: $(test_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) $(CXXFLAGS) *.cpp > .depend

-include .depend
