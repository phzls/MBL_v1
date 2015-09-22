//
// Created by Liangsheng Zhang on 9/17/15.
//

#ifndef MBL_V1_HAM_EVOL_MODEL_H
#define MBL_V1_HAM_EVOL_MODEL_H

/*
 * This file contains particular models of Hamiltonian operators
 */

#include <sstream>
#include <iostream>
#include <vector>
#include "ham_evol.h"

using namespace std;

/*
 * This operator constructs the Heisenberg random Hamiltonian with random fields
 * given in cosine form. It only constructs one sector with given total spin Z
 * The total spin Z is passed in as a parameter, where up spin is considered 1
 * and down spin is considered -1
 */
class HamEvolHeisenRandomCosSzSector : public HamEvolVanillaReal
{
private:
    const double h_; // Disorder strength

    const int total_spin_z_; // Total spin z value, where up spin is 1 and down spin is -1

    // Spin Configuration. It gives the number in full basic basis for each index in
    // our restricted basic basis. The numbers will be sorted in ascending order
    vector<int> config_;

    void Spin_Config_();
    bool config_constructed_; // Whether spin configuration has been constructed

    void Repr_Init_(); // Initialize the representation string stream as well as type

    void Dim_Cal_(); // Calculate the true dimension of the restricted Hilbert space

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field

public:
    HamEvolHeisenRandomCosSzSector(int size, double h, int total_spin_z = 0, bool debug = false):
            HamEvolVanillaReal(size), h_(h), total_spin_z_(total_spin_z), debug_(debug),
            config_constructed_(false) { Repr_Init_(); Dim_Cal_(); }

    // Construct the Hamiltonian only for the sector of specifiec total Z spin,
    // in the basis of restricted basic states whose total Z spin satisfy the condition.
    // The map between the index in the restricted basis and full basic basis is given
    // in config_
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Construct Transition Matrix
    void Transition_Compute(TransitionMatrix&, const string&) const {
        cout << "Not implemented yet." << endl;
        abort();
    }

    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXcd>&) const {
        cout << "Not implemented yet." << endl;
        abort();
    }

    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXd>&) const {
        cout << "Not implemented yet." << endl;
        abort();
    }

    // Return config_ so that one can map from restricted binary basis to full binary basis
    // This is introduced as a hack as it doesn't fit into Transition_Compute paragon
    vector<int> Get_Spin_Config() const {return config_;}

    string Eigen_Basis_Type() const {return "Restricted Basic";}

    virtual ~HamEvolHeisenRandomCosSzSector() {};
};








// ====================================================================================================








/*
 * This operator constructs the Heisenberg Hamiltonian with quasi-periodic fields
 * given in cosine form. It only constructs one sector with given total spin Z
 * The total spin Z is passed in as a parameter, where up spin is considered 1
 * and down spin is considered -1
 */
class HamEvolHeisenQuasiSzSector : public HamEvolVanillaReal
{
private:
    const double h_; // Disorder strength

    const int total_spin_z_; // Total spin z value, where up spin is 1 and down spin is -1

    // Spin Configuration. It gives the number in full basic basis for each index in
    // our restricted basic basis. The numbers will be sorted in ascending order
    vector<int> config_;

    void Spin_Config_();
    bool config_constructed_; // Whether spin configuration has been constructed

    void Repr_Init_(); // Initialize the representation string stream as well as type

    void Dim_Cal_(); // Calculate the true dimension of the restricted Hilbert space

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field

public:
    HamEvolHeisenQuasiSzSector(int size, double h, int total_spin_z = 0, bool debug = false):
            HamEvolVanillaReal(size), h_(h), total_spin_z_(total_spin_z), debug_(debug),
            config_constructed_(false) { Repr_Init_(); Dim_Cal_(); }

    // Construct the Hamiltonian only for the sector of specifiec total Z spin,
    // in the basis of restricted basic states whose total Z spin satisfy the condition.
    // The map between the index in the restricted basis and full basic basis is given
    // in config_
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Construct Transition Matrix
    void Transition_Compute(TransitionMatrix&, const string&) const {
        cout << "Not implemented yet." << endl;
        abort();
    }

    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXcd>&) const {
        cout << "Not implemented yet." << endl;
        abort();
    }

    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXd>&) const {
        cout << "Not implemented yet." << endl;
        abort();
    }

    // Return config_ so that one can map from restricted binary basis to full binary basis
    // This is introduced as a hack as it doesn't fit into Transition_Compute paragon
    vector<int> Get_Spin_Config() const {return config_;}

    string Eigen_Basis_Type() const {return "Restricted Basic";}

    virtual ~HamEvolHeisenQuasiSzSector() {};
};

#endif //MBL_V1_HAM_EVOL_MODEL_H
