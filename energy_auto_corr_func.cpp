//
// Created by Liangsheng Zhang on 9/25/15.
//

/*
 * This file implements some functions related to energy_auto_corr
 */

#include "energy_auto_corr_func.h"

using namespace std;
using namespace Eigen;

void Ham_for_Ising_Random_Simp_Shift_Cos_Real_Flo(const EvolOP* floquet, MatrixXcd& ham){
    floquet -> Get_Ham(ham, "Eigen", "");
}

void Ham_for_Heisen_Random_Cos_Sz_Sector_Shift_Real_Flo(const EvolOP* floquet, MatrixXcd& ham){
    floquet -> Get_Ham(ham, "Eigen", "");
}

