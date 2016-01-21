//
// Created by Liangsheng Zhang on 9/28/15.
//

/**
 ** Implementation of FloEvolHeisenQuasiSzSectorShiftRealTau class
 **/

#include <complex>
#include <iostream>
#include <algorithm>
#include <map>
#include "constants.h"
#include "flo_evol_model.h"
#include "screen_output.h"
#include "randomc.h"
#include "combinatorics.h"

using namespace std;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolHeisenQuasiSzSectorShiftRealTau::Repr_Init_(){
    repr_ << "Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Floquet_L=" << size_ << ",h=" << h_ << ",tau=" << tau_;
    type_ = "Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Floquet";
}

/*
 * Initialize quasi-periodic fields
 */
void FloEvolHeisenQuasiSzSectorShiftRealTau::Evol_Para_Init() {
    random_h_.resize(size_);
    double u = RanGen_mersenne.Random();
    double theta = 2 * Pi * u; // Random initial phase

    // Random Fields
    for (int i=0; i<size_;i++){
        random_h_[i] = h_ * cos(i * Pi / Phi + theta);
    }

    if (debug_){
        cout << "Random phase: " << theta << endl;
        cout << "Longitude fields:" << endl;
        for (int i=0; i<size_; i++) cout << random_h_[i] << endl;
        cout << endl;
    }
}

/*
 * Copy quasi-periodic fields
 */
void FloEvolHeisenQuasiSzSectorShiftRealTau::Evol_Para_Copy(const vector< vector<double> >& random_h) {
    random_h_.resize(size_);

    cout << "Evol_Para_Copy not tested!!" << endl;
    abort();

    if(random_h.size() != 1){
        cout << type_ << " has only one set of quasi-periodic fields." << endl;
        cout << "Obtained: " << random_h.size() << endl;
        abort();
    }

    if(random_h[0].size() != size_){
        cout << "Incorrect number of quasi-periodic fields passed in." << endl;
        cout << "Expected: " << size_ << " Obtained: " << random_h[0].size() << endl;
        abort();
    }

    // Copy quasi-periodic fields and check the bounds
    const double lower_bound = -h_;
    const double upper_bound = h_;
    for(int i=0; i<size_; i++){
        double val = random_h[0][i];
        if(val < lower_bound || val > upper_bound){
            cout << "Quasi-periodic field passed in at pos " << i << " is outside of bound." << endl;
            cout << "Lower bound: " << lower_bound << " Upper bound: " << upper_bound << endl;
            cout << "Val: " << val << endl;
        }
        random_h_[i] = val;
    }

    if (debug_){
        cout << "Quasi-periodic fields:" << endl;
        for (int i=0; i<size_;i++) cout << random_h_[i] << endl;
    }
}

/*
 * Compute the dimension of the restricted Hilbert space for that sector
 */
void FloEvolHeisenQuasiSzSectorShiftRealTau::Dim_Cal_() {
    int temp = total_spin_z_ + size_;
    if( temp % 2 != 0){
        cout << "Cannot achieve total Z spin " << total_spin_z_ << " for system size " << size_ << endl;
        abort();
    }

    int up_spin = temp/2; // Number of up spins
    dim_ = binomial_coef(size_, up_spin);
}

/*
 * Construct availabe spin configurations: it gives basic basis numbers whose total
 * Z spin agrees with the specified number
 * The config_ vector will be kept even when the Hamiltonian is destroyed to used
 * for possible future basis transformation
 */

void FloEvolHeisenQuasiSzSectorShiftRealTau::Spin_Config_() {
    if(!config_constructed_) { // Otherwise no need to construct
        int temp = total_spin_z_ + size_;
        if( temp % 2 != 0){
            cout << "Cannot achieve total Z spin " << total_spin_z_ << " for system size " << size_ << endl;
            abort();
        }

        int up_spin = temp/2; // Number of up spins
        comb_num(0, size_-1, up_spin, config_);

        if(config_.size() != dim_){
            cout << "Incompatible dimension for " << repr_ << endl;
            cout << "Number of configurations: " << config_.size() << endl;
            cout << "Dimension of restricted Hilbert space: " << dim_ << endl;
            abort();
        }

        config_constructed_ = true;

        // Sort in ascending order
        sort(config_.begin(), config_.end());

        if(debug_){
            cout << "Spin states in basic basis with total z spin " << total_spin_z_ << ":" << endl;
            for(int i=0; i<config_.size(); i++)
                cout << config_[i] << endl;
        }
    }
}

/*
 * Construct the time evolution matrix
 */
void FloEvolHeisenQuasiSzSectorShiftRealTau::Evol_Construct() {

    if (!constructed_){
        evol_op_ = MatrixXcd::Zero(dim_, dim_);
        evol_op_real_ = MatrixXd::Zero(dim_, dim_);
        evol_op_imag_ = MatrixXd::Zero(dim_, dim_);
        constructed_ = true;
    }
    else{
        cout << "Evolution operator has been constructed." << endl;
        abort();
    }

    if(size_ < 1) {
        cout << "Cannot construct Floquet with " << size_ << " number of spins for " << type_ << endl;
        abort();
    }

    Spin_Config_();

    MatrixXcd evol_half_x = MatrixXcd::Zero(dim_,dim_);
    MatrixXcd evol_z = MatrixXcd::Zero(dim_,dim_);

    Evol_X_Z_Construct_(evol_half_x, evol_z);

    evol_op_ = evol_half_x * evol_z * evol_half_x;

    for (int i=0; i<evol_op_.cols();i++){
        for (int j=0; j<evol_op_.rows();j++){
            evol_op_real_(j,i) = real(evol_op_(j,i));
            evol_op_imag_(j,i) = imag(evol_op_(j,i));
        }
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<evol_op_real_.cols(); i++){
        for (int j=i+1; j< evol_op_real_.rows();j++){
            if (abs(evol_op_real_(j,i) - evol_op_real_(i,j)) > 5.0e-6){
                cout << "Real part is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << evol_op_real_(i,j) << endl;
                cout << "(j,i): " << evol_op_real_(j,i) << endl;
                abort();
            }
        }
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<evol_op_imag_.cols(); i++){
        for (int j=i+1; j< evol_op_imag_.rows();j++){
            if (abs(evol_op_imag_(j,i) - evol_op_imag_(i,j)) > 5.0e-6){
                cout << "Imaginary part is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << evol_op_imag_(i,j) << endl;
                cout << "(j,i): " << evol_op_imag_(j,i) << endl;
                abort();
            }
        }
    }

    if (debug_){
        cout << "Sqrt X part of evolution matrix:" << endl;
        matrix_write(evol_half_x);
        cout << endl;

        cout << "Z part of evolution matrix:" << endl;
        matrix_write(evol_z);
        cout << endl;

        cout << "Evolution matrix:" << endl;
        matrix_write(evol_op_);
        cout << endl;

        cout << "Real part of evolution matrix:" << endl;
        matrix_write(evol_op_real_);
        cout << endl;

        cout << "Imaginary part of evolution matrix:" << endl;
        matrix_write(evol_op_imag_);
        cout << endl;
    }
}

/*
 * Construct the Hamiltonian only for the sector of given total Z spin
 */
void FloEvolHeisenQuasiSzSectorShiftRealTau::Evol_X_Z_Construct_(MatrixXcd& evol_half_x, MatrixXcd& evol_z) {
    // Change config_ to a map for easy lookup of basic states to restricted states
    map<int, int> config_map;
    for(int i=0; i<config_.size(); i++) {
        int state = config_[i];
        config_map[state] = i;
    }

    for(int i=0; i<config_.size(); i++){
        int state = config_[i];
        int mod_state = state; // The state that would be modified later
        int spin = 0;
        int prev_spin = 0;
        int flip_state = -1; // The state after some spins are flipped in full basic basis
        int restricted_flip_state = -1; // Corresponding state in restricted basic basis
        // for flip_state

        int pos = 1; // Representing the position of spin
        int prev_pos = 1; // Representing the position of prev_spin

        double z_val = 0; // diagonal term

        // The leftmost spin
        spin = 2*(mod_state & 1) - 1;
        mod_state = mod_state >> 1;
        z_val += random_h_[0] * spin;

        for(int j=1; j<size_; j++){
            prev_spin = spin;
            prev_pos = pos;

            spin = 2*(mod_state & 1) - 1;
            mod_state = mod_state >> 1;
            pos = pos << 1;

            z_val += random_h_[j] * spin;
            z_val += spin * prev_spin;

            if(prev_spin + spin == 0){
                // We should stay in the same sector of given total Z spin
                // so the flipped two spins must point to opposite directions
                flip_state = state ^ prev_pos; // flip prev_spin
                flip_state = flip_state ^ pos; // flip spin

                restricted_flip_state = config_map[flip_state];
                evol_half_x(restricted_flip_state, i) += 1;
            }
            evol_z(i, i) = exp( complex<double>(0,-tau_*z_val) );
        }
    }

    // Check if the x_half matrix is Hermitian
    for (int i=0; i<evol_half_x.cols(); i++){
        for (int j=i+1; j< evol_half_x.rows();j++){
            if (abs(evol_half_x(j,i) - evol_half_x(i,j)) > 1.0e-10){
                cout << "half_x Hamiltonian is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << evol_half_x(i,j) << endl;
                cout << "(j,i): " << evol_half_x(j,i) << endl;
                abort();
            }
        }
    }

    SelfAdjointEigenSolver<MatrixXcd> x_eigen;

    x_eigen.compute(evol_half_x);

    for (int i=0; i<dim_; i++){
        for (int j=0; j< dim_; j++){
            if (i==j){
                evol_half_x(i,j) = exp( - complex<double>(tau_,0) * Complex_I * x_eigen.eigenvalues()[i] );
            }
            else{
                evol_half_x(i,j) = complex<double>(0,0);
            }
        }
    }

    evol_half_x = x_eigen.eigenvectors() * evol_half_x * x_eigen.eigenvectors().adjoint();
}

/*
 * Construct the following Hamiltonian regardless of what input s is:
 * sum_i { sigma_x^i*sigma_x^(i+1) + sigma_y^i*sigma_y^(i+1) + sigma_z^i*sigma_z^(i+1) + h_i*sigma_x^i }
 */
void FloEvolHeisenQuasiSzSectorShiftRealTau::Get_Ham(MatrixXcd& ham, string basis, string s) const {
    ham = MatrixXcd::Zero(dim_, dim_);

    if(!config_constructed_){
        cout << "Spin configuration is not initialized. Cannot construct Hamiltonian directly" << endl;
        abort();
    }

    // Change config_ to a map for easy lookup of basic states to restricted states
    map<int, int> config_map;
    for(int i=0; i<config_.size(); i++) {
        int state = config_[i];
        config_map[state] = i;
    }

    for(int i=0; i<config_.size(); i++){
        int state = config_[i];
        int mod_state = state; // The state that would be modified later
        int spin = 0;
        int prev_spin = 0;
        int flip_state = -1; // The state after some spins are flipped in full basic basis
        int restricted_flip_state = -1; // Corresponding state in restricted basic basis
        // for flip_state

        int pos = 1; // Representing the position of spin
        int prev_pos = 1; // Representing the position of prev_spin

        double z_val = 0; // diagonal term

        // The leftmost spin
        spin = 2*(mod_state & 1) - 1;
        mod_state = mod_state >> 1;
        z_val += random_h_[0] * spin;

        for(int j=1; j<size_; j++){
            prev_spin = spin;
            prev_pos = pos;

            spin = 2*(mod_state & 1) - 1;
            mod_state = mod_state >> 1;
            pos = pos << 1;

            z_val += random_h_[j] * spin;
            z_val += spin * prev_spin;

            if(prev_spin + spin == 0){
                // We should stay in the same sector of given total Z spin
                // so the flipped two spins must point to opposite directions
                flip_state = state ^ prev_pos; // flip prev_spin
                flip_state = flip_state ^ pos; // flip spin

                restricted_flip_state = config_map[flip_state];
                ham(restricted_flip_state, i) += 2;
            }
            ham(i, i) = z_val;
        }
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<ham.cols(); i++){
        for (int j=i+1; j< ham.rows();j++){
            if (abs(ham(j,i) - ham(i,j)) > 1.0e-10){
                cout << "Hamiltonian is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << ham(i,j) << endl;
                cout << "(j,i): " << ham(j,i) << endl;
                abort();
            }
        }
    }

    if (debug_){
        cout << "Hamiltonian in restricted basic basis:" << endl;
        complex_matrix_write(ham);
        cout << endl;
    }

    if(basis == "Restricted Basic") return;
    else if (basis == "Eigen"){
        // Change to eigenvector basis
        if(!eigen_info_){
            cout << "Hamiltonian under eigen basis for " << repr_
            << " cannot be established if eigenvectors are not computed." << endl;
            abort();
        }

        if(debug_){
            cout << "Eigenvector:" << endl;
            real_matrix_write(evec_);
            cout << endl;
        }

        ham = evec_.adjoint() * ham * evec_;
    }
    else{
        cout << "Hamiltonion under " << basis << " cannot be established for " << repr_ << endl;
        abort();
    }

    if(debug_){
        cout << "Hamiltonain under " << basis << " basis:" << endl;
        complex_matrix_write(ham);
        cout << endl;
    }
}




