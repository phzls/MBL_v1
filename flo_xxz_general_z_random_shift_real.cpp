//
// Created by Liangsheng Zhang on 12/4/15.
//

/*
 * Implementation of FloEvolXXZGeneralZRandomShiftReal class
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
void FloEvolXXZGeneralZRandomShiftReal::Repr_Init_(){
    // Check the rv_type
    if(rv_type_ != "Gaussian" && rv_type_ != "Uniform"){
        cout << rv_type_ << " random variables cannot be used" << endl;
        abort();
    }

    repr_ << "XXZ_" << rv_type_ << "_Random_Shift_Real_Floquet_L=" << size_
          << ",sigma=" << sigma_;
    type_ = "XXZ_" + rv_type_ + "_Random_Shift_Real_Floquet";
}

/*
 * Initialize random numbers
 */
void FloEvolXXZGeneralZRandomShiftReal::Evol_Para_Init() {
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

    if (debug_){
        cout << rv_type_ << " random numbers:" << endl;
        for (int i=0; i<size_;i++) cout << random_[i] << endl;
    }

    initialized_ = true;
}

/*
 * Construct the time evolution matrix
 */
void FloEvolXXZGeneralZRandomShiftReal::Evol_Construct() {
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
        cout << "Gaussians for " << repr_ << " have not been constructed." << endl;
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
void FloEvolXXZGeneralZRandomShiftReal::Evol_X_Construct_(MatrixXcd & evol_half_x) {
    vector<double>cos_sin_prod(size_+1); // Index means the number of different spins

    const double factor = sqrt(1-sigma_*sigma_);

    cos_sin_prod[0] = 1;

    // Compute the sine part
    for (int i=1; i<size_+1; i++){
        cos_sin_prod[i] = cos_sin_prod[i-1] * ( -sin( (tau_*g_*factor)/2.0 ) );
    }

    // Add the cosine part
    double cos_prod = 1;
    for (int i=size_; i>=0; i-- ){
        cos_sin_prod[i] *= cos_prod;
        cos_prod *= cos( (g_*tau_*factor)/2.0 );
    }

    for (int i=0; i<dim_; i++){
        for (int j=i; j<dim_;j++){
            int flip = __builtin_popcount(i ^ j);
            complex<double> sign;
            if (flip % 4 == 1) sign = Complex_I;
            else if (flip % 4 == 2) sign = -Complex_one;
            else if (flip % 4 == 3) sign = -Complex_I;
            else sign = Complex_one;

            evol_half_x(j,i) = sign * complex<double>(cos_sin_prod.at(flip), 0);
            if (i != j) evol_half_x(i,j) = evol_half_x(j,i);
        }
    }
}

/*
 *  Construct z part of the evolution matrix
 */
void FloEvolXXZGeneralZRandomShiftReal::Evol_Z_Construct_(MatrixXcd & evol_z) {
    const int int_size = sizeof(int)*CHAR_BIT; // Number of bits in an integer
    for (int i=0; i<dim_;i++){
        bitset<int_size> state(i);
        double value = 0;
        for (int j=0; j<size_;j++){
            value += (h_+g_*sigma_*random_[j])*(2*state[j]-1);
            if(j>0) value += (2*state[j]-1)*(2*state[j-1]-1);
        }

        evol_z(i,i) = exp( complex<double>(0,-tau_*value) );
    }
}

/*
 * Construct the following Hamiltonian no matter what the input string is
 */
void FloEvolXXZGeneralZRandomShiftReal::Get_Ham(MatrixXcd& ham, string basis, string s){
    ham = MatrixXcd::Zero(dim_, dim_);
    const int int_size = sizeof(int)*CHAR_BIT; // Number of bits in an integer
    const double factor = sqrt(1-sigma_*sigma_);
    // Construct the matrix in basic basis
    for(int i=0; i<dim_; i++){
        bitset<int_size> state(i);

        double z_value = 0; // diagonal element

        for (int j=0; j<size_;j++){
            z_value += (h_+g_*sigma_*random_[j])*(2*state[j]-1);
            if(j>0) z_value += (2*state[j]-1)*(2*state[j-1]-1);

            state.flip(j);
            int x_state = state.to_ulong();
            state.flip(j); // Flip back the spin

            ham(x_state, i) += g_*factor;
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



