//
// Created by Liangsheng Zhang on 4/15/15.
//

#include <iostream>
#include <string>
#include "model_func.h"
#include "flo_evol_model.h"

using namespace std;

// For Ising random simple floquet operator
string Flo_Evol_Ising_Random_Simp_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingRandomSimp(size, W, debug);

    string type = model -> Type();
    replace(type.begin(), type.end(), '_', ' ');

    return type;
}

// For Ising random simple shift real floquet operator
string Flo_Evol_Ising_Random_Simp_Shift_Real_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingRandomSimpShiftReal(size, W, debug);

    string type = model -> Type();
    replace(type.begin(), type.end(), '_', ' ');

    return type;
}