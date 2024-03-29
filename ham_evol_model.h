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










// =============================================================================================================











/*
 * This operator constructs the Ising simple random Hamiltonian corresponding to
 * FloEvolIsingRandomSimpShiftCosReal
 */
class HamEvolIsingRandomSimpCos : public HamEvolVanillaReal
{
private:
    const double W_; // Disorder strength

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field

public:
    HamEvolIsingRandomSimpCos(int size, double W, bool debug = false): HamEvolVanillaReal(size), W_(W), debug_(debug)
    { Repr_Init_();}

    // Construct the Hamiltonian
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    string Eigen_Basis_Type() const {return "Basic";}

    virtual ~HamEvolIsingRandomSimpCos() {};
};









// ===========================================================================================================









/*
 * This operator constructs the Ising simple quasi-periodic Hamiltonian corresponding to
 * FloEvolIsingQuasiSimpShiftCosReal
 */
class HamEvolIsingQuasiSimp : public HamEvolVanillaReal
{
private:
    const double W_; // Disorder strength

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field

public:
    HamEvolIsingQuasiSimp(int size, double W, bool debug = false): HamEvolVanillaReal(size), W_(W), debug_(debug)
    { Repr_Init_();}

    // Construct the Hamiltonian
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    string Eigen_Basis_Type() const {return "Basic";}

    virtual ~HamEvolIsingQuasiSimp() {};
};








// ==============================================================================================================







/*
 * This operator constructs the modified Heisenberg random Hamiltonian with random
 * fields given in cosine form. It only constructs one sector with given total spin
 * Z The total spin Z is passed in as a parameter, where up spin is considered 1
 * and down spin is considered -1
 */
class HamEvolHeisenModifiedRandomCosSzSector : public HamEvolVanillaReal
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
    HamEvolHeisenModifiedRandomCosSzSector(int size, double h, int total_spin_z = 0, bool debug = false):
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

    virtual ~HamEvolHeisenModifiedRandomCosSzSector() {};
};









// =======================================================================================================









/*
 * This operator constructs the modified Heisenberg Hamiltonian with quasi-periodic fields
 * given in cosine form. It only constructs one sector with given total spin Z
 * The total spin Z is passed in as a parameter, where up spin is considered 1
 * and down spin is considered -1
 */
class HamEvolHeisenModifiedQuasiSzSector : public HamEvolVanillaReal
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
    HamEvolHeisenModifiedQuasiSzSector(int size, double h, int total_spin_z = 0, bool debug = false):
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

    virtual ~HamEvolHeisenModifiedQuasiSzSector() {};
};









// ==============================================================================================================







/*
 * This operator constructs the continuous modified Heisenberg random Hamiltonian
 * with random fields given in cosine form. The "continuous" refers to (1-alpha*h)
 * discount of zz term. xx and yy terms are discounted by 1-h
 * This Hamiltonian only constructs one sector with given total spin Z The total
 * spin Z is passed in as a parameter, where up spin is considered 1 and down spin
 * is considered -1
 */
class HamEvolHeisenConModifiedRandomCosSzSector : public HamEvolVanillaReal
{
private:
    const double h_; // Disorder strength

    const double alpha_; // The parameter for zz discount

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
    HamEvolHeisenConModifiedRandomCosSzSector(int size, double h, double alpha, int total_spin_z = 0,
                                              bool debug = false):
            HamEvolVanillaReal(size), h_(h), alpha_(alpha), total_spin_z_(total_spin_z),
            debug_(debug), config_constructed_(false) { Repr_Init_(); Dim_Cal_(); }

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

    virtual ~HamEvolHeisenConModifiedRandomCosSzSector() {};
};









// =======================================================================================================









/*
 * This operator constructs the modified Heisenberg Hamiltonian with quasi-periodic fields
 * given in cosine form. The "continuous" refers to (1-alpha*h) discount of zz term.
 * xx and yy terms are discounted by 1-h
 * It only constructs one sector with given total spin Z The total spin Z is passed in as
 * a parameter, where up spin is considered 1 and down spin is considered -1
 */
class HamEvolHeisenConModifiedQuasiSzSector : public HamEvolVanillaReal
{
private:
    const double h_; // Disorder strength

    const double alpha_; // The parameter for zz discount

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
    HamEvolHeisenConModifiedQuasiSzSector(int size, double h, double alpha,
                                          int total_spin_z = 0, bool debug = false):
            HamEvolVanillaReal(size), h_(h), alpha_(alpha), total_spin_z_(total_spin_z),
            debug_(debug), config_constructed_(false) { Repr_Init_(); Dim_Cal_(); }

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

    virtual ~HamEvolHeisenConModifiedQuasiSzSector() {};
};



#endif //MBL_V1_HAM_EVOL_MODEL_H
