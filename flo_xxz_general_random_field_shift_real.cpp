//
// Created by Liangsheng Zhang on 12/7/15.
//

/*
 * Implementation of FloEvolXXZGeneralRandomFieldShiftReal class
 */

#include <complex>
#include <iostream>
#include <cmath>
#include <bitset>
#include <climits>
#include "constants.h"
#include "flo_evol_model.h"
#include "screen_output.h"
#include "randomc.h"

using namespace std;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolXXZGeneralRandomFieldShiftReal::Repr_Init_(){
    // Check the rv_type
    if(rv_type_ != "Gaussian" && rv_type_ != "Uniform"){
        cout << rv_type_ << " random variables cannot be used" << endl;
        abort();
    }

    repr_ << "XXZ_" << rv_type_ << "_Random_Field_Shift_Real_Floquet_L=" << size_
    << ",sigma=" << sigma_;
    type_ = "XXZ_" + rv_type_ + "_Random_Field_Shift_Real_Floquet";
}

/*
 * Initialize random numbers
 */
void FloEvolXXZGeneralRandomFieldShiftReal::Evol_Para_Init() {
    random_.resize(size_);

    // Random Fields
    if(rv_type_ == "Gaussian"){
        for (int i=0; i<size_;){
            // Using Box Muller transformation to general standard normal
            double u1 = RanGen_mersenne.Random();
            double u2 = RanGen_mersenne.Random();
            random_[i] = sqrt(-2*log(u1))*cos(2*Pi*u2);
            i++;
            if(i<size_) random_[i] = sqrt(-2*log(u1))*sin(2*Pi*u2);
            i++;
        }
    }
    else if(rv_type_ == "Uniform"){
        for (int i=0; i<size_; i++){
            double u = RanGen_mersenne.Random();
            random_[i] = 2*sqrt(3)*u - sqrt(3);
        }
    }

    J_ = sqrt( 1 - sigma_*sigma_ * ( (h_*h_ + g_*g_)/( 1 + g_*g_ + h_*h_ ) ) );

    if (debug_){
        cout << rv_type_ << " random numbers:" << endl;
        for (int i=0; i<size_;i++) cout << random_[i] << endl;
        cout << "sigma: " << sigma_ << endl;
        cout << "J: " << J_ << endl;
    }

    initialized_ = true;
}

/*
 * Copy random numbers
 */
void FloEvolXXZGeneralRandomFieldShiftReal::Evol_Para_Copy(const vector< vector<double> >& random) {
    random_.resize(size_);

    cout << "Evol_Para_Copy not tested!!" << endl;
    abort();

    if(random.size() != 1){
        cout << type_ << " has only one set of random fields." << endl;
        cout << "Obtained: " << random.size() << endl;
        abort();
    }

    if(random[0].size() != size_){
        cout << "Incorrect number of random fields passed in." << endl;
        cout << "Expected: " << size_ << " Obtained: " << random[0].size() << endl;
        abort();
    }

    // Copy Random Fields and check the bounds if there are
    if(rv_type_ == "Gaussian"){
        random_ = random[0];
    }
    else if(rv_type_ == "Uniform"){
        const double lower_bound = -sqrt(3);
        const double upper_bound = sqrt(3);
        for(int i=0; i<size_; i++){
            double val = random[0][i];
            if(val < lower_bound || val > upper_bound){
                cout << "Random field passed in at pos " << i << " is outside of bound." << endl;
                cout << "Lower bound: " << lower_bound << " Upper bound: " << upper_bound << endl;
                cout << "Val: " << val << endl;
                abort();
            }
            random_[i] = val;
        }
    }

    if (debug_){
        cout << "sigma: " << sigma_ << endl;
        cout << rv_type_ << " random numbers:" << endl;
        for (int i=0; i<size_;i++) cout << random_[i] << endl;
    }

    initialized_ = true;
}

/*
 * Construct the time evolution matrix
 */
void FloEvolXXZGeneralRandomFieldShiftReal::Evol_Construct() {
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

    if(!initialized_){
        cout << rv_type_ << " for " << repr_ << " have not been constructed." << endl;
        abort();
    }

    MatrixXcd evol_half_x = MatrixXcd::Zero(dim_,dim_);
    MatrixXcd evol_z = MatrixXcd::Zero(dim_,dim_);

    Evol_X_Construct_(evol_half_x);
    Evol_Z_Construct_(evol_z);

    evol_op_ = evol_half_x * evol_z * evol_half_x;

    evol_op_real_ = evol_op_.real();
    evol_op_imag_ = evol_op_.imag();

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
 * Construct X part of the evolution matrix
 */
void FloEvolXXZGeneralRandomFieldShiftReal::Evol_X_Construct_(MatrixXcd & evol_half_x) {
    const double sq = sqrt(h_*h_+g_*g_);
    for(int i=0; i<dim_; i++){

        for(int j=0; j<dim_; j++){
            complex<double> prod = 1;
            int state1 = i;
            int state2 = j;

            for(int k=0; k<size_; k++){
                int spin1 = 2*(state1 & 1) - 1;
                int spin2 = 2*(state2 & 1) - 1;
                state1 = state1 >> 1;
                state2 = state2 >> 1;

                double factor = (J_ + sigma_*random_[k]) * sq * tau_ / 2.0;

                if(spin1 == spin2) prod *= complex<double>( cos(factor),
                                                            - h_ / sq * spin1 * sin(factor) );
                else prod *= complex<double>( 0, -g_/sq*sin(factor) );
            }

            evol_half_x(i,j) = prod;
            if(i!=j) evol_half_x(j,i) = prod;
        }
    }
}

/*
 *  Construct z part of the evolution matrix
 */
void FloEvolXXZGeneralRandomFieldShiftReal::Evol_Z_Construct_(MatrixXcd & evol_z) {
    for (int i=0; i<dim_;i++){
        double value = 0;
        int state = i;
        int prev_spin = 2*(state & 1) - 1;

        for (int j=1; j<size_; j++){
            state = state >> 1;
            int spin = 2*(state & 1) - 1;
            value += J_*spin*prev_spin;
            prev_spin = spin;
        }

        evol_z(i,i) = exp( complex<double>(0,-tau_*value) );
    }
}

/*
 * Construct the following Hamiltonian no matter what the input string is
 */
void FloEvolXXZGeneralRandomFieldShiftReal::Get_Ham(MatrixXcd& ham, string basis, string s){
    ham = MatrixXcd::Zero(dim_, dim_);
    const int int_size = sizeof(int)*CHAR_BIT; // Number of bits in an integer

    // Construct the matrix in basic basis
    for(int i=0; i<dim_; i++){
        bitset<int_size> state(i);

        double z_value = 0; // diagonal element

        for (int j=0; j<size_;j++){
            z_value += h_*(J_+sigma_*random_[j])*(2*state[j]-1);
            if(j>0) z_value += J_*(2*state[j]-1)*(2*state[j-1]-1);

            state.flip(j);
            int x_state = state.to_ulong();
            state.flip(j); // Flip back the spin

            ham(x_state, i) += g_*(J_+sigma_*random_[j]);
        }

        ham(i, i) = z_value;
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<ham.cols(); i++){
        for (int j=i+1; j< ham.rows();j++){
            if (abs(ham(j,i) - ham(i,j)) > 5.0e-6){
                cout << "Hamiltonian is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << ham(i,j) << endl;
                cout << "(j,i): " << ham(j,i) << endl;
                abort();
            }
        }
    }

    if(basis == "Basic") return;
    else if (basis == "Eigen"){
        // Change to eigenvector basis
        if(!eigen_info_){
            cout << "Hamiltonian under eigen basis for " << repr_
            << " cannot be established if eigenvectors are not computed." << endl;
            abort();
        }

        ham = evec_.adjoint() * ham * evec_;
    }
    else{
        cout << "Hamiltonion under " << basis << " cannot be established for " << repr_ << endl;
        abort();
    }
}





