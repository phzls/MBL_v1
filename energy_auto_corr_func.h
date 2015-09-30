//
// Created by Liangsheng Zhang on 9/25/15.
//

#ifndef MBL_V1_ENERGY_AUTO_CORR_FUNC_H
#define MBL_V1_ENERGY_AUTO_CORR_FUNC_H

#include <Eigen/Dense>
#include "evol_op.h"

void Ham_for_Ising_Random_Simp_Shift_Cos_Real_Flo(const EvolOP*, Eigen::MatrixXcd&);

void Ham_for_Heisen_Random_Cos_Sz_Sector_Shift_Real_Flo(const EvolOP*, MatrixXcd&);

void Ham_for_Heisen_Quasi_Sz_Sector_Shift_Real_Flo(const EvolOP*, MatrixXcd&);

void Ham_for_Heisen_Random_Cos_Sz_Sector_Modified_Flo(const EvolOP*, MatrixXcd&);

#endif //MBL_V1_ENERGY_AUTO_CORR_FUNC_H
