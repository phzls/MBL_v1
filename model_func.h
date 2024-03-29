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

// For Heisen random cos sz sector Hamiltonian operator
string Ham_Evol_Heisen_Random_Cos_Sz_Sector_Func(const AllPara&, EvolOP*&);

// For Heisen quasi-periodic sz sector Hamiltonian operator
string Ham_Evol_Heisen_Quasi_Sz_Sector_Func(const AllPara&, EvolOP*&);

// For Ising random simple shift cos real tau floquet operator
string Flo_Evol_Ising_Random_Simp_Shift_Cos_Real_Tau_Func(const AllPara&, EvolOP*&);

// For Heisen random cos sz sector real shift tau floquet operator
string Flo_Evol_Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Func(const AllPara&, EvolOP*&);

// For Heisen quasi-periodic sz sector real shift tau floquet operator
string Flo_Evol_Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Func(const AllPara&, EvolOP*&);

// For Ising random simple cos hamiltonian operator
string Ham_Evol_Ising_Random_Simp_Cos_Func(const AllPara&, EvolOP*&);

// For Ising quasi-periodic simple hamiltonian operator
string Ham_Evol_Ising_Quasi_Simp_Func(const AllPara&, EvolOP*&);

// For Heisen random cos sz sector modified tau floquet operator
string Flo_Evol_Heisen_Random_Cos_Sz_Sector_Modified_Tau_Func(const AllPara&, EvolOP*&);

// For Heisen quasi-periodic sz sector modified tau floquet operator
string Flo_Evol_Heisen_Quasi_Sz_Sector_Modified_Tau_Func(const AllPara&, EvolOP*&);

// For Heisen modified random cos sz sector Hamiltonian operator
string Ham_Evol_Heisen_Modified_Random_Cos_Sz_Sector_Func(const AllPara&, EvolOP*&);

// For Heisen modified quasi-periodic sz sector Hamiltonian operator
string Ham_Evol_Heisen_Modified_Quasi_Sz_Sector_Func(const AllPara&, EvolOP*&);

// For Heisen continuous modified random cos sz sector Hamiltonian operator
string Ham_Evol_Heisen_Con_Modified_Random_Cos_Sz_Sector_Func(const AllPara&, EvolOP*&);

// For Heisen continuous modified quasi-periodic sz sector Hamiltonian operator
string Ham_Evol_Heisen_Con_Modified_Quasi_Sz_Sector_Func(const AllPara&, EvolOP*&);

// For XXZ gaussian random shift real floquet operator with Gaussian random fields
string Flo_Evol_XXZ_Gaussian_Random_Shift_Real_Func(const AllPara&, EvolOP*&);

// For XXZ uniform random shift real floquet operator with uniform random fields
string Flo_Evol_XXZ_Uniform_Random_Shift_Real_Func(const AllPara&, EvolOP*&);

// For XXZ gaussian random Z shift real floquet operator with Gaussian random fields
// along z direction
string Flo_Evol_XXZ_Gaussian_Z_Random_Shift_Real_Func(const AllPara&, EvolOP*&);

// For XXZ uniform random Z shift real floquet operator with uniform random fields
// along z direction
string Flo_Evol_XXZ_Uniform_Z_Random_Shift_Real_Func(const AllPara&, EvolOP*&);

// For XXZ gaussian random field shift real floquet operator with Gaussian random fields
// along z direction
string Flo_Evol_XXZ_Gaussian_Random_Field_Shift_Real_Func(const AllPara&, EvolOP*&);

// For XXZ uniform random field shift real floquet operator with uniform random fields
// along z direction
string Flo_Evol_XXZ_Uniform_Random_Field_Shift_Real_Func(const AllPara&, EvolOP*&);

#endif //MBL_V1_MODEL_FUNC_H
