//
// Created by Liangsheng Zhang on 4/14/15.
//

#ifndef MBL_V1_FLO_EVOL_MODEL_H
#define MBL_V1_FLO_EVOL_MODEL_H

/*
 * This file contains particular models of Floquet time evolution operators
 */

#include <sstream>
#include <iostream>
#include <vector>
#include "constants.h"
#include "flo_evol.h"

using namespace std;

/*
 * This operator constructs the flo_evol_xxz_random_simp which is a simplified version of
 * flo_evol_xxz_random with many parameters equal to 1
 */
class FloEvolIsingRandomSimp : public FloEvolVanilla
{
private:
    const double W_; // Disorder strength

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field part at each site

public:
    FloEvolIsingRandomSimp(int size, double W, bool debug = false):
            FloEvolVanilla(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // No parameter to initialize
    void Evol_Para_Init();

    virtual ~FloEvolIsingRandomSimp() {};
};







//=======================================================================================================






/*
 * This operator constructs the shifted version of flo_evol_xxz_random_simp which is a simplified version of
 * flo_evol_xxz_random with many parameters equal to 1. Effectively, it is sqrt(U_x) * U_z * sqrt(U_x). It
 * utilizes the feature that this operator has real eigenvectors, so only real eigenvectors are returned
 * and only real part of eigenvalues are computed.
 */
class FloEvolIsingRandomSimpShiftReal : public FloEvolVanillaReal
{
private:
    const double W_; // Disorder strength

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field part at each site

public:
    FloEvolIsingRandomSimpShiftReal(int size, double W, bool debug = false):
            FloEvolVanillaReal(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // No parameter to initialize
    void Evol_Para_Init();

    virtual ~FloEvolIsingRandomSimpShiftReal() {};
};









//=======================================================================================================






/*
 * This operator constructs the shifted version of flo_evol_xxz_random_simp but using a quasi-periodic
 * longitude field instead of a truly random field
 */
class FloEvolIsingQuasiSimpShiftReal : public FloEvolVanillaReal
{
private:
    const double W_; // Disorder strength

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field part at each site

public:
    FloEvolIsingQuasiSimpShiftReal(int size, double W, bool debug = false):
            FloEvolVanillaReal(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // No parameter to initialize
    void Evol_Para_Init();

    virtual ~FloEvolIsingQuasiSimpShiftReal() {};
};


#endif //MBL_V1_FLO_EVOL_MODEL_H
