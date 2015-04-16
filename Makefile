INC = ~/

models = flo_evol.o flo_ising_random_simp.o flo_ising_random_simp_shift_real.o

tasks = disorder_transition.o

disorder_transition = zz_corr_square.o disorder_model_transition.o zz_time_corr.o zz_all_time_corr.o zz_time_corr_component.o ent_var.o

small_methods = matrix_algebra.o screen_output.o mersenne.o evec_to_basic.o reduced_density_left_2.o generic_func.o

transitions = basis_transition.o basic_full.o basic_parity.o parity_full.o

test_objects = test.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions)

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