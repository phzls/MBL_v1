//
// Created by Liangsheng Zhang on 4/15/15.
//

#include <iostream>
#include <string>
#include "model_func.h"
#include "flo_evol_model.h"
#include "ham_evol_model.h"

using namespace std;

// For Ising random simple floquet operator
string Flo_Evol_Ising_Random_Simp_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingRandomSimp(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Ising random simple shift real floquet operator
string Flo_Evol_Ising_Random_Simp_Shift_Real_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingRandomSimpShiftReal(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Ising quasi-periodic simple shift real floquet operator
string Flo_Evol_Ising_Quasi_Simp_Shift_Real_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingQuasiSimpShiftReal(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Ising all random simple shift real floquet operator
string Flo_Evol_Ising_All_Random_Simp_Shift_Real_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingAllRandomSimpShiftReal(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Ising all quasi simple shift real floquet operator
string Flo_Evol_Ising_All_Quasi_Simp_Shift_Real_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingAllQuasiSimpShiftReal(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Ising random simple shift real floquet operator
string Flo_Evol_Ising_Random_Simp_Shift_Cos_Real_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingRandomSimpShiftCosReal(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Heisen random cos sz sector Hamiltonian operator
string Ham_Evol_Heisen_Random_Cos_Sz_Sector_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // random field strength
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new HamEvolHeisenRandomCosSzSector(size, h, total_spin_z, debug);

    string type = model -> Type();

    return type;
}


// For Heisen quasi-periodic sz sector Hamiltonian operator
string Ham_Evol_Heisen_Quasi_Sz_Sector_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // random field strength
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new HamEvolHeisenQuasiSzSector(size, h, total_spin_z, debug);

    string type = model -> Type();

    return type;
}

// For Ising random simple shift cos real tau floquet operator
string Flo_Evol_Ising_Random_Simp_Shift_Cos_Real_Tau_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength
    const double tau = parameters.floquet.tau; // Period

    const bool debug = parameters.generic.debug;

    model = new FloEvolIsingRandomSimpShiftCosRealTau(size, W, tau, debug);

    string type = model -> Type();

    return type;
}

// For Heisen random cos sz sector real shift tau floquet operator
string Flo_Evol_Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // Disorder strength
    const double tau = parameters.floquet.tau; // Period
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new FloEvolHeisenRandomCosSzSectorShiftRealTau(size, h, tau, total_spin_z, debug);

    string type = model -> Type();

    return type;
}

// For Heisen quasi-periodic sz sector real shift tau floquet operator
string Flo_Evol_Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // Disorder strength
    const double tau = parameters.floquet.tau; // Period
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new FloEvolHeisenQuasiSzSectorShiftRealTau(size, h, tau, total_spin_z, debug);

    string type = model -> Type();

    return type;
}

// For Ising random simple cos hamiltonian operator
string Ham_Evol_Ising_Random_Simp_Cos_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new HamEvolIsingRandomSimpCos(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Ising quasi-periodic simple hamiltonian operator
string Ham_Evol_Ising_Quasi_Simp_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double W = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new HamEvolIsingQuasiSimp(size, W, debug);

    string type = model -> Type();

    return type;
}

// For Heisen random cos sz sector modified tau floquet operator
string Flo_Evol_Heisen_Random_Cos_Sz_Sector_Modified_Tau_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // Disorder strength
    const double tau = parameters.floquet.tau; // Period
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new FloEvolHeisenRandomCosSzSectorModifiedTau(size, h, tau, total_spin_z, debug);

    string type = model -> Type();

    return type;
}

// For Heisen quasi-periodic sz sector modified tau floquet operator
string Flo_Evol_Heisen_Quasi_Sz_Sector_Modified_Tau_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // Disorder strength
    const double tau = parameters.floquet.tau; // Period
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new FloEvolHeisenQuasiSzSectorModifiedTau(size, h, tau, total_spin_z, debug);

    string type = model -> Type();

    return type;
}

// For Heisen modified random cos sz sector Hamiltonian operator
string Ham_Evol_Heisen_Modified_Random_Cos_Sz_Sector_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // random field strength
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new HamEvolHeisenModifiedRandomCosSzSector(size, h, total_spin_z, debug);

    string type = model -> Type();

    return type;
}


// For Heisen modified quasi-periodic sz sector Hamiltonian operator
string Ham_Evol_Heisen_Modified_Quasi_Sz_Sector_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // random field strength
    const int total_spin_z = parameters.floquet.total_spin_z;

    const bool debug = parameters.generic.debug;

    model = new HamEvolHeisenModifiedQuasiSzSector(size, h, total_spin_z, debug);

    string type = model -> Type();

    return type;
}

// For Heisen constinuous modified random cos sz sector Hamiltonian operator
string Ham_Evol_Heisen_Con_Modified_Random_Cos_Sz_Sector_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // random field strength
    const int total_spin_z = parameters.floquet.total_spin_z;
    const double alpha = parameters.floquet.alpha;

    const bool debug = parameters.generic.debug;

    model = new HamEvolHeisenConModifiedRandomCosSzSector(size, h, alpha, total_spin_z, debug);

    string type = model -> Type();

    return type;
}


// For Heisen continuous modified quasi-periodic sz sector Hamiltonian operator
string Ham_Evol_Heisen_Con_Modified_Quasi_Sz_Sector_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double h = parameters.floquet.J; // random field strength
    const int total_spin_z = parameters.floquet.total_spin_z;
    const double alpha = parameters.floquet.alpha;

    const bool debug = parameters.generic.debug;

    model = new HamEvolHeisenConModifiedQuasiSzSector(size, h, alpha, total_spin_z, debug);

    string type = model -> Type();

    return type;
}


// For XXZ gaussian random shift real floquet operator with Gaussian random fields
string Flo_Evol_XXZ_Gaussian_Random_Shift_Real_Func(const AllPara& parameters, EvolOP*& model){
    const int size = parameters.generic.size; // System Size
    const double J = parameters.floquet.J; // Disorder strength

    const bool debug = parameters.generic.debug;

    model = new FloEvolXXZGaussianRandomShiftReal(size, J, debug);

    string type = model -> Type();

    return type;
}

