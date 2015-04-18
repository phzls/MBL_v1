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
void Ising_Random_Simp_Flo_Para(AllPara&);

// For Ising random simple shift real floquet operator
void Ising_Random_Simp_Shift_Real_Flo_Para(AllPara&);

// For Ising quasi-periodic simple shift real floquet operator
void Ising_Quasi_Simp_Shift_Real_Flo_Para(AllPara&);

// For Ising all random simple shift real floquet operator
void Ising_All_Random_Simp_Shift_Real_Flo_Para(AllPara&);


//=================================== TASKS =====================================================

// generic parameters
void Generic_Para(AllPara&);

// Output parameters
void Output_Para(AllPara&);

// For disorder_transition
void Disorder_Transition_Para(AllPara&);



#endif //MBL_V1_PARA_MODEL_TASK_H
