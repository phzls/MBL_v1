INC = ~/

models = flo_evol.o flo_ising_random_simp.o flo_ising_random_simp_shift_real.o flo_ising_quasi_simp_shift_real.o flo_ising_all_random_simp_shift_real.o flo_ising_all_quasi_simp_shift_real.o flo_ising_random_simp_shift_cos_real.o ham_evol.o ham_evol_heisen_random_cos_sz_sector.o ham_evol_heisen_quasi_sz_sector.o flo_ising_random_simp_shift_cos_real_tau.o flo_heisen_random_cos_sz_sector_shift_real_tau.o flo_heisen_quasi_sz_sector_shift_real_tau.o ham_evol_ising_random_simp_cos.o ham_evol_ising_quasi_simp.o flo_heisen_random_cos_sz_sector_modified_tau.o flo_heisen_quasi_sz_sector_modified_tau.o ham_evol_heisen_modified_quasi_sz_sector.o ham_evol_heisen_modified_random_cos_sz_sector.o

tasks = disorder_transition.o single_model_time_evolution.o state_evol.o multi_model_time_evolution.o density_evol.o

disorder_transition = zz_corr_square.o disorder_model_transition.o zz_time_corr.o zz_all_corr_square.o zz_time_corr_component.o ent_var.o zz_corr_square_all_sample.o zz_time_corr_all_sample.o log_zz_all_corr_square.o log_zz_corr_square.o zz_all_time_corr.o ent_var_all_mean.o flo_level_stats_ave.o ent_scaled_mean.o ham_level_stats_ave.o

model_evolution = init_obj.o evol_data.o random_pure_init.o product_random_init.o random_product_init.o entropy_per_model.o init_obj_func.o left_spin_random_init.o leftmost_spin_per_model.o

op_auto_corr = flo_op_auto_corr.o op_auto_corr_func.o op_corr_local_para.o energy_auto_corr.o energy_auto_corr_func.o

controls = tasks_models.o model_func.o

small_methods = matrix_algebra.o screen_output.o mersenne.o evec_to_basic.o reduced_density_left_2.o generic_func.o combinatorics.o disorder_transition_func.o

transitions = basis_transition.o basic_full.o basic_parity.o parity_full.o

para_get = para_read.o para_model_task.o

test_objects = test.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls) $(model_evolution) $(op_auto_corr)
mbl_test_objects = main_test.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls) $(model_evolution) $(op_auto_corr)
mbl_objects = main.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls) $(model_evolution) $(op_auto_corr)
task_model_objects = task_type_print.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls) $(model_evolution) $(op_auto_corr)
auto_objects = main_auto.o $(models) $(tasks) $(disorder_transition) $(small_methods) $(transitions) $(controls) $(para_get) $(model_evolution) $(op_auto_corr)

CXXFLAGS = -O3 -fopenmp
CXX = g++

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
CXX = icpc
CXXFLAGS += -I $(INC)
endif

OUT = flo_xxz_all_random_simp_shift_real_12

all: mbl_real

auto: $(auto_objects)
	$(CXX) -o $(OUT) $(CXXFLAGS) $^

mbl_real: $(mbl_objects)
	$(CXX) -o $(OUT) $(CXXFLAGS) $^

mbl: $(mbl_test_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

test: $(test_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

task_model: $(task_model_objects)
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) $(CXXFLAGS) *.cpp > .depend

-include .depend
