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

    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >&);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

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

    vector<double> random_h_; // Random longitude field

public:
    FloEvolIsingRandomSimpShiftReal(int size, double W, bool debug = false):
            FloEvolVanillaReal(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >&);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string s = "");

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

    vector<double> random_h_; // Quasi-Random longitude field

public:
    FloEvolIsingQuasiSimpShiftReal(int size, double W, bool debug = false):
            FloEvolVanillaReal(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize quasi-periodic fields
    void Evol_Para_Init();

    // Copy quasi-periodic fields
    void Evol_Para_Copy(const vector< vector<double> >&);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real quasi-periodic fields to copy." << endl;
        abort();
    }

    virtual ~FloEvolIsingQuasiSimpShiftReal() {};
};










//=======================================================================================================









/*
 * This operator constructs the shifted version of flo_evol_xxz_random_simp where all interaction
 * strengths are random
 */
class FloEvolIsingAllRandomSimpShiftReal : public FloEvolVanillaReal
{
private:
    const double W_; // Disorder strength

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field
    vector<double> random_J_; // Random nearest neighbor interaction
    vector<double> random_g_; // Random transverse field

public:
    FloEvolIsingAllRandomSimpShiftReal(int size, double W, bool debug = false):
            FloEvolVanillaReal(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    virtual ~FloEvolIsingAllRandomSimpShiftReal() {};
};











//=======================================================================================================









/*
 * This operator constructs the shifted version of flo_evol_xxz_random_simp where all interaction
 * strengths are quasi-periodic
 */
class FloEvolIsingAllQuasiSimpShiftReal : public FloEvolVanillaReal
{
private:
    const double W_; // Disorder strength

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field
    vector<double> random_J_; // Random nearest neighbor interaction
    vector<double> random_g_; // Random transverse field

public:
    FloEvolIsingAllQuasiSimpShiftReal(int size, double W, bool debug = false):
            FloEvolVanillaReal(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize quasi-periodic fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    virtual ~FloEvolIsingAllQuasiSimpShiftReal() {};
};










//=============================================================================================================










/*
 * This operator constructs the shifted version of flo_evol_xxz_random_simp but using a cosine with
 * phases as the longitude field
 */
class FloEvolIsingRandomSimpShiftCosReal : public FloEvolVanillaReal
{
private:
    const double W_; // Disorder strength

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field

public:
    FloEvolIsingRandomSimpShiftCosReal(int size, double W, bool debug = false):
            FloEvolVanillaReal(size), W_(W), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    virtual ~FloEvolIsingRandomSimpShiftCosReal() {};
};








// ========================================================================================





/*
 * This operator constructs the shifted version of flo_evol_xxz_random_simp using a cosine with
 * phases as the longitude field. The period can also be varied
 */
class FloEvolIsingRandomSimpShiftCosRealTau : public FloEvolVanillaReal
{
private:
    const double W_; // Disorder strength
    const double tau_; // Period

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_h_; // Random longitude field

public:
    FloEvolIsingRandomSimpShiftCosRealTau(int size, double W, double tau, bool debug = false):
            FloEvolVanillaReal(size), W_(W), tau_(tau), debug_(debug) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string) const;

    virtual ~FloEvolIsingRandomSimpShiftCosRealTau() {};
};











// ===========================================================================================================











/*
 * This operator constructs the floquet corresponding to HamEvolHeisenRandomCosSzSector
 * Particularly, since Sz is still conserved, only one sector is constructed
 * The operator is still shifted to make sure that eigenvectors are real
 */
class FloEvolHeisenRandomCosSzSectorShiftRealTau : public FloEvolVanillaReal
{
private:
    const double h_; // Disorder strength
    const double tau_; // Period
    const int total_spin_z_; // Total spin z value, where up spin is 1 and down spin is -1

    // Construct x and z part of time evolution operator
    void Evol_X_Z_Construct_(MatrixXcd&, MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    void Dim_Cal_(); // Calculate the true dimension of the restricted Hilbert space

    const bool debug_; // Used for debug outputs

    void Spin_Config_();
    bool config_constructed_; // Whether spin configuration has been constructed

    vector<double> random_h_; // Random longitude field

    // Spin Configuration. It gives the number in full basic basis for each index in
    // our restricted basic basis. The numbers will be sorted in ascending order
    vector<int> config_;

public:
    FloEvolHeisenRandomCosSzSectorShiftRealTau(int size, double h, double tau, int total_spin_z = 0,
                                               bool debug = false):  FloEvolVanillaReal(size), h_(h), tau_(tau),
                                                                     total_spin_z_(total_spin_z), debug_(debug),
                                                                     config_constructed_(false)
    { Repr_Init_(); Dim_Cal_(); }

    // Construct evolutionary operator only for the sector of specifiec total Z spin,
    // in the basis of restricted basic states whose total Z spin satisfy the condition.
    // The map between the index in the restricted basis and full basic basis is given
    // in config_
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

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

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string) const;

    // Return config_ so that one can map from restricted binary basis to full binary basis
    // This is introduced as a hack as it doesn't fit into Transition_Compute paragon
    vector<int> Get_Spin_Config() const {return config_;}

    string Eigen_Basis_Type() const {return "Restricted Basic";}

    virtual ~FloEvolHeisenRandomCosSzSectorShiftRealTau() {};
};







// =========================================================================================================








/*
 * This operator constructs the floquet corresponding to HamEvolHeisenQuasiSzSector
 * Particularly, since Sz is still conserved, only one sector is constructed
 * The operator is still shifted to make sure that eigenvectors are real
 */
class FloEvolHeisenQuasiSzSectorShiftRealTau : public FloEvolVanillaReal
{
private:
    const double h_; // Disorder strength
    const double tau_; // Period
    const int total_spin_z_; // Total spin z value, where up spin is 1 and down spin is -1

    // Construct x and z part of time evolution operator
    void Evol_X_Z_Construct_(MatrixXcd&, MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    void Dim_Cal_(); // Calculate the true dimension of the restricted Hilbert space

    const bool debug_; // Used for debug outputs

    void Spin_Config_();
    bool config_constructed_; // Whether spin configuration has been constructed

    vector<double> random_h_; // Random longitude field

    // Spin Configuration. It gives the number in full basic basis for each index in
    // our restricted basic basis. The numbers will be sorted in ascending order
    vector<int> config_;

public:
    FloEvolHeisenQuasiSzSectorShiftRealTau(int size, double h, double tau, int total_spin_z = 0,
                                               bool debug = false):  FloEvolVanillaReal(size), h_(h), tau_(tau),
                                                                     total_spin_z_(total_spin_z), debug_(debug),
                                                                     config_constructed_(false)
    { Repr_Init_(); Dim_Cal_(); }

    // Construct evolutionary operator only for the sector of specifiec total Z spin,
    // in the basis of restricted basic states whose total Z spin satisfy the condition.
    // The map between the index in the restricted basis and full basic basis is given
    // in config_
    void Evol_Construct();

    // Initialize quasi-periodic fields
    void Evol_Para_Init();

    // Copy quasi-periodic fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

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

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string) const;

    // Return config_ so that one can map from restricted binary basis to full binary basis
    // This is introduced as a hack as it doesn't fit into Transition_Compute paragon
    vector<int> Get_Spin_Config() const {return config_;}

    string Eigen_Basis_Type() const {return "Restricted Basic";}

    virtual ~FloEvolHeisenQuasiSzSectorShiftRealTau() {};
};













// ===========================================================================================================











/*
 * A modified version of FloEvolHeisenRandomCosSzSectorShiftRealTau. The xx and yy terms
 * are discounted by (1-h) so that the system can localize when h approaches 1
 */
class FloEvolHeisenRandomCosSzSectorModifiedTau : public FloEvolVanillaReal
{
private:
    const double h_; // Disorder strength
    const double tau_; // Period
    const int total_spin_z_; // Total spin z value, where up spin is 1 and down spin is -1

    // Construct x and z part of time evolution operator
    void Evol_X_Z_Construct_(MatrixXcd&, MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    void Dim_Cal_(); // Calculate the true dimension of the restricted Hilbert space

    const bool debug_; // Used for debug outputs

    void Spin_Config_();
    bool config_constructed_; // Whether spin configuration has been constructed

    vector<double> random_h_; // Random longitude field

    // Spin Configuration. It gives the number in full basic basis for each index in
    // our restricted basic basis. The numbers will be sorted in ascending order
    vector<int> config_;

public:
    FloEvolHeisenRandomCosSzSectorModifiedTau(int size, double h, double tau, int total_spin_z = 0,
                                               bool debug = false):  FloEvolVanillaReal(size), h_(h), tau_(tau),
                                                                     total_spin_z_(total_spin_z), debug_(debug),
                                                                     config_constructed_(false)
    { Repr_Init_(); Dim_Cal_(); }

    // Construct evolutionary operator only for the sector of specifiec total Z spin,
    // in the basis of restricted basic states whose total Z spin satisfy the condition.
    // The map between the index in the restricted basis and full basic basis is given
    // in config_
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

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

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string) const;

    // Return config_ so that one can map from restricted binary basis to full binary basis
    // This is introduced as a hack as it doesn't fit into Transition_Compute paragon
    vector<int> Get_Spin_Config() const {return config_;}

    string Eigen_Basis_Type() const {return "Restricted Basic";}

    virtual ~FloEvolHeisenRandomCosSzSectorModifiedTau() {};
};







// =========================================================================================================








/*
 * A modified version of FloEvolHeisenQuasiCosSzSectorShiftRealTau. The xx and yy terms
 * are discounted by (1-h) so that the system can localize when h approaches 1
 */
class FloEvolHeisenQuasiSzSectorModifiedTau : public FloEvolVanillaReal
{
private:
    const double h_; // Disorder strength
    const double tau_; // Period
    const int total_spin_z_; // Total spin z value, where up spin is 1 and down spin is -1

    // Construct x and z part of time evolution operator
    void Evol_X_Z_Construct_(MatrixXcd&, MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    void Dim_Cal_(); // Calculate the true dimension of the restricted Hilbert space

    const bool debug_; // Used for debug outputs

    void Spin_Config_();
    bool config_constructed_; // Whether spin configuration has been constructed

    vector<double> random_h_; // Random longitude field

    // Spin Configuration. It gives the number in full basic basis for each index in
    // our restricted basic basis. The numbers will be sorted in ascending order
    vector<int> config_;

public:
    FloEvolHeisenQuasiSzSectorModifiedTau(int size, double h, double tau, int total_spin_z = 0,
                                           bool debug = false):  FloEvolVanillaReal(size), h_(h), tau_(tau),
                                                                 total_spin_z_(total_spin_z), debug_(debug),
                                                                 config_constructed_(false)
    { Repr_Init_(); Dim_Cal_(); }

    // Construct evolutionary operator only for the sector of specifiec total Z spin,
    // in the basis of restricted basic states whose total Z spin satisfy the condition.
    // The map between the index in the restricted basis and full basic basis is given
    // in config_
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

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

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string) const;

    // Return config_ so that one can map from restricted binary basis to full binary basis
    // This is introduced as a hack as it doesn't fit into Transition_Compute paragon
    vector<int> Get_Spin_Config() const {return config_;}

    string Eigen_Basis_Type() const {return "Restricted Basic";}

    virtual ~FloEvolHeisenQuasiSzSectorModifiedTau() {};
};









// ========================================================================================================










/*
 * This operator constructs the shifted version of flo_evol_xxz_random which is full version of
 * flo_evol_xxz_random with Gaussian random fields along particular direction. The trace of
 * the corresponding H^2 is also fixed, which determines the relation between sigma and J.
 *
 * Effectively, it is sqrt(U_x) * U_z * sqrt(U_x). It utilizes the feature that this operator
 * has real eigenvectors, so only real eigenvectors are returned and only real part of
 * eigenvalues are computed.
 *
 * Currently the model is parametrized by sigma
 */
class FloEvolXXZGaussianRandomShiftReal : public FloEvolVanillaReal
{
private:
    const double h_; // Longitude field
    const double g_; // Transverse field
    const double tau_; // Time step

    double J_; // Disorder strength
    double sigma_; // The SD for Gaussian random variable

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> gaussian_; // A set of standard normal random variables

    bool initialized_; // Whether the gaussians and sigma_ have been computed

public:
    /*FloEvolXXZGaussianRandomShiftReal(int size, double J, bool debug = false):
            FloEvolVanillaReal(size), J_(J), debug_(debug), h_(0.8090),
            g_(0.9045), tau_(0.8) { Repr_Init_();}*/

    FloEvolXXZGaussianRandomShiftReal(int size, double sigma, bool debug = false):
            FloEvolVanillaReal(size), sigma_(sigma), debug_(debug), h_(0.8090),
            g_(0.9045), tau_(0.8) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string s = "");

    virtual ~FloEvolXXZGaussianRandomShiftReal() {};
};












// =============================================================================================











/*
 * This operator constructs the shifted version of flo_evol_xxz_random which is full version of
 * flo_evol_xxz_random with uniform random fields along particular direction. The trace of
 * the corresponding H^2 is also fixed, which determines the relation between sigma and J.
 *
 * Effectively, it is sqrt(U_x) * U_z * sqrt(U_x). It utilizes the feature that this operator
 * has real eigenvectors, so only real eigenvectors are returned and only real part of
 * eigenvalues are computed.
 *
 * Currently the model is parametrized by sigma
 */
class FloEvolXXZUniformRandomShiftReal : public FloEvolVanillaReal
{
private:
    const double h_; // Longitude field
    const double g_; // Transverse field
    const double tau_; // Time step

    double J_; // Disorder strength
    double sigma_; // The SD for uniform random variable

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> uniform_; // A set of uniform random variables on -sqrt(3) to sqrt(3)

    bool initialized_; // Whether the gaussians and sigma_ have been computed

public:
    /*FloEvolXXZUniformRandomShiftReal(int size, double J, bool debug = false):
            FloEvolVanillaReal(size), J_(J), debug_(debug), h_(0.8090),
            g_(0.9045), tau_(0.8) { Repr_Init_();}*/

    FloEvolXXZUniformRandomShiftReal(int size, double sigma, bool debug = false):
            FloEvolVanillaReal(size), sigma_(sigma), debug_(debug), h_(0.8090),
            g_(0.9045), tau_(0.8) { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string s = "");

    virtual ~FloEvolXXZUniformRandomShiftReal() {};
};









// ===============================================================================================







/*
 * This operator constructs the shifted version of flo_evol_xxz_random which is full version of
 * flo_evol_xxz_random with generic random fields along z direction. The trace of the corresponding
 * H^2 is fixed. The random variables all have mean 0 and variance 1. Its type is passed in
 * as arguments during construction
 *
 * Effectively, it is sqrt(U_x) * U_z * sqrt(U_x). It utilizes the feature that this operator
 * has real eigenvectors, so only real eigenvectors are returned and only real part of
 * eigenvalues are computed.
 *
 * Currently the model is parametrized by sqrt(1-sigma^2)
 */
class FloEvolXXZGeneralZRandomShiftReal : public FloEvolVanillaReal
{
private:
    const double h_; // Longitude field
    const double g_; // Transverse field
    const double tau_; // Time step

    const double sigma_; // The SD for random variable

    const string rv_type_; // The type of random variable

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_; // A set of random variables with mean 0 and variance 1

    bool initialized_; // Whether the random variables and sigma_ have been computed

public:
    FloEvolXXZGeneralZRandomShiftReal(int size, double transverse, string rv_type = "Gaussian",
                                      bool debug = false): FloEvolVanillaReal(size),
                                                           sigma_(sqrt(1-transverse*transverse)),
                                                           debug_(debug), h_(0.8090), g_(0.9045),
                                                           tau_(0.8), rv_type_(rv_type)
    { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string s = "");

    virtual ~FloEvolXXZGeneralZRandomShiftReal() {};
};









// ===============================================================================================







/*
 * This operator constructs the shifted version of flo_evol_xxz_random_field which is full
 * version of flo_evol_xxz_random with generic random fields along (g,0,h) direction.
 * The trace of the corresponding H^2 is fixed in the thermodynamic limit, which gives a
 * relation between J and sigma, the SD of random fields. The random variables constructed
 * all have mean 0 and variance 1. Its type is passed in as arguments during construction
 *
 * Effectively, it is sqrt(U_f) * U_z * sqrt(U_f), where U_f is the exponential of the whole
 * random field part of the Hamiltonian. It utilizes the feature that this operator
 * has real eigenvectors, so only real eigenvectors are returned and only real part of
 * eigenvalues are computed.
 *
 * Currently the model is parametrized by sigma
 */
class FloEvolXXZGeneralRandomFieldShiftReal : public FloEvolVanillaReal
{
private:
    const double h_; // Longitude field
    const double g_; // Transverse field
    const double tau_; // Time step

    double sigma_; // The SD for random variable
    double J_; // The strenght of zz coupling

    const string rv_type_; // The type of random variable

    // Construct x part of time evolution operator
    void Evol_X_Construct_(MatrixXcd&);

    // Construct z part of time evolution operator
    void Evol_Z_Construct_(MatrixXcd&);

    void Repr_Init_(); // Initialize the representation string stream as well as type

    const bool debug_; // Used for debug outputs

    vector<double> random_; // A set of random variables with mean 0 and variance 1

    bool initialized_; // Whether the random variables and sigma_ have been computed

public:
    FloEvolXXZGeneralRandomFieldShiftReal(int size, double sigma, string rv_type = "Gaussian",
                                      bool debug = false): FloEvolVanillaReal(size), sigma_(sigma),
                                                           debug_(debug), h_(0.8090), g_(0.9045),
                                                           tau_(0.8), rv_type_(rv_type)
    { Repr_Init_();}

    // Construct evolutionary operator
    void Evol_Construct();

    // Initialize random fields
    void Evol_Para_Init();

    // Copy random fields
    void Evol_Para_Copy(const vector< vector<double> >& para);
    void Evol_Para_Copy(const vector< vector< complex<double> > >& para){
        cout << type_ << " must have real random fields to copy." << endl;
        abort();
    }

    // Construct the Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    void Get_Ham(MatrixXcd&, string, string s = "");

    virtual ~FloEvolXXZGeneralRandomFieldShiftReal() {};
};




#endif //MBL_V1_FLO_EVOL_MODEL_H
