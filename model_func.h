//
// Created by Liangsheng Zhang on 4/15/15.
//

/*
 * This file contains functions which generating a pointer to a model. Type of the model is
 * returned as a string
 */

#ifndef MBL_V1_MODEL_FUNC_H
#define MBL_V1_MODEL_FUNC_H

#include "evol_op.h"
#include "parameters.h"

using namespace std;

// For Ising random simple floquet operator
string Flo_Evol_Ising_Random_Simp_Func(const AllPara&, EvolOP*&);

// For Ising random simple shift real floquet operator
string Flo_Evol_Ising_Random_Simp_Shift_Real_Func(const AllPara&, EvolOP*&);

// For Ising quasi-periodic simple shift real floquet operator
string Flo_Evol_Ising_Quasi_Simp_Shift_Real_Func(const AllPara&, EvolOP*&);

// For Ising all random simple shift real floquet operator
string Flo_Evol_Ising_All_Random_Simp_Shift_Real_Func(const AllPara&, EvolOP*&);

// For Ising all quasi-periodic simple shift real floquet operator
string Flo_Evol_Ising_All_Quasi_Simp_Shift_Real_Func(const AllPara&, EvolOP*&);

// For Ising random simple shift cos real floquet operator
string Flo_Evol_Ising_Random_Simp_Shift_Cos_Real_Func(const AllPara&, EvolOP*&);

#endif //MBL_V1_MODEL_FUNC_H
