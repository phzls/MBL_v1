//
// Created by Liangsheng Zhang on 4/18/15.
//

/*
 * Functions related to reading parameters from files for all models and tasks
 */

#ifndef MBL_V1_PARA_MODEL_TASK_H
#define MBL_V1_PARA_MODEL_TASK_H

#include "para_read.h"

using namespace std;


//=================================== MODELS ===================================================

// For Ising random simple floquet operator
void Ising_Random_Simp_Flo_Para(AllPara&, string);

// For Ising random simple shift real floquet operator
void Ising_Random_Simp_Shift_Real_Flo_Para(AllPara&, string);

// For Ising quasi-periodic simple shift real floquet operator
void Ising_Quasi_Simp_Shift_Real_Flo_Para(AllPara&, string);

// For Ising all random simple shift real floquet operator
void Ising_All_Random_Simp_Shift_Real_Flo_Para(AllPara&, string);

// For Ising all quasi-periodic simple shift real floquet operator
void Ising_All_Quasi_Simp_Shift_Real_Flo_Para(AllPara&, string);

// For Ising random simple shift cosine real floquet operator
void Ising_Random_Simp_Shift_Cos_Real_Flo_Para(AllPara&, string);

// For Heisenberg random cosine Sz sector Hamiltonian Operator
void Heisen_Random_Cos_Sz_Sector_Ham_Para(AllPara&, string);

// For Heisenberg quasi-periodic cosine Sz sector Hamiltonian Operator
void Heisen_Quasi_Sz_Sector_Ham_Para(AllPara&, string);

// For Ising random simple shift cosine real tau floquet operator
void Ising_Random_Simp_Shift_Cos_Real_Tau_Flo_Para(AllPara&, string);

// For Heisenberg random cos sz sector shift real tau floquet operator
void Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Flo_Para(AllPara&, string);

// For Heisenberg quasi-periodic sz sector shift real tau floquet operator
void Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Flo_Para(AllPara&, string);

// For Ising random simple cosine hamiltonian operator
void Ising_Random_Simp_Cos_Ham_Para(AllPara&, string);

// For Ising quasi-periodic simple hamiltonian operator
void Ising_Quasi_Simp_Ham_Para(AllPara&, string);


//=================================== TASKS =====================================================

// generic parameters
void Generic_Para(AllPara&, string);

// Output parameters
void Output_Para(AllPara&, string);

// For disorder_transition
void Disorder_Transition_Para(AllPara&, string);

// For single_model_time_evolution_para
void Single_Model_Time_Evolution_Para(AllPara&, string);

// For multi_model_time_evolution_para
void Multi_Model_Time_Evolution_Para(AllPara&, string);

// For operator autocorrelation. Now only working with Floquet dynamics
void Op_Auto_Corr_Para(AllPara&, string);

#endif //MBL_V1_PARA_MODEL_TASK_H
