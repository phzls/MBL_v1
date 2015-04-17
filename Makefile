INC = ~/

models = flo_evol.o flo_ising_random_simp.o flo_ising_random_simp_shift_real.o flo_ising_quasi_simp_shift_real.o

tasks = disorder_transition.o

disorder_transition = zz_corr_square.o disorder_model_transition.o zz_time_corr.o zz_all_time_corr.o zz_time_corr_component.o ent_var.o

controls = tasks_models.o model_func.o

small_methods = matrix_algebra.o screen_output.o mersenne.o evec_to_basic.o reduced_density_left_2.o generic_func.o

transitions = basis_transition.o basic_full.o basic_parity.o parity_full.o

test_objects = test.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls)
mbl_test_objects = main_test.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls)
mbl_objects = main.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls)

CXXFLAGS = -O3 -fopenmp
CXX = g++

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
CXX = icpc
CXXFLAGS += -I $(INC)
endif

OUT = flo_xxz_random_simp_shift_real_11

all: mbl_real

mbl_real: $(mbl_objects)
	$(CXX) -o $(OUT) $(CXXFLAGS) $^

mbl: $(mbl_test_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

test: $(test_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) $(CXXFLAGS) *.cpp > .depend

-include .depend
