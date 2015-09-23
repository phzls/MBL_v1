//
// Created by Liangsheng Zhang on 9/23/15.
//

/*
 * Implementation of FloEvolIsingQuasiSimpCosRealTau class
 */

#include <complex>
#include <iostream>
#include <cmath>
#include "constants.h"
#include "flo_evol_model.h"
#include "screen_output.h"
#include "randomc.h"

using namespace std;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolIsingRandomSimpShiftCosRealTau::Repr_Init_(){
    repr_ << "Ising_Random_Simp_Shift_Cos_Real_Tau_Floquet_L=" << size_ << ",W=" << W_ << ",tau=" << tau_;
    type_ = "Ising_Random_Simp_Shift_Cos_Real_Tau_Floquet";
}

/*
 * Initialize random numbers
 */
void FloEvolIsingRandomSimpShiftCosRealTau::Evol_Para_Init() {
    random_h_.resize(size_);
    vector<double> theta(size_);

    // Random Fields
    for (int i=0; i<size_;i++){
        // Random phase
        double u = RanGen_mersenne.Random();
        theta[i] = 2 * Pi * u;
        random_h_[i] = 1 + W_ * cos(theta[i]);
    }

    if (debug_){
        cout << "Random phase: " << endl;
        for (int i=0; i<size_; i++) cout << theta[i] << endl;
        cout << "Longitude fields:" << endl;
        for (int i=0; i<size_; i++) cout << random_h_[i] << endl;
        cout << endl;
    }
}

/*
 * Construct the time evolution matrix
 */
void FloEvolIsingRandomSimpShiftCosRealTau::Evol_Construct() {

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

    MatrixXcd evol_half_x = MatrixXcd::Zero(dim_,dim_);
    MatrixXcd evol_z = MatrixXcd::Zero(dim_,dim_);

    Evol_X_Construct_(evol_half_x);
    Evol_Z_Construct_(evol_z);

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
 * Construct X part of the evolution matrix
 */
void FloEvolIsingRandomSimpShiftCosRealTau::Evol_X_Construct_(MatrixXcd & evol_half_x) {
    vector<double>cos_prod(dim_+1);
    vector<double>sin_prod(dim_+1);

    cos_prod[0] = 1;
    sin_prod[0] = 1;

    for (int i=1; i<dim_+1; i++){
        cos_prod[i] = cos_prod[i-1] * cos((1-W_)/2);
        sin_prod[i] = sin_prod[i-1] * ( - sin((1-W_)/2) );
    }

    for (int i=0; i<dim_; i++){
        for (int j=i; j<dim_;j++){
            int flip = __builtin_popcount(i ^ j);
            complex<double> sign;
            if (flip % 4 == 1) sign = Complex_I;
            else if (flip % 4 == 2) sign = -Complex_one;
            else if (flip % 4 == 3) sign = -Complex_I;
            else sign = Complex_one;

            evol_half_x(j,i) = sign * complex<double>(tau_ * cos_prod[size_ - flip] * sin_prod[flip], 0);
            if (i != j) evol_half_x(i,j) = evol_half_x(j,i);
        }
    }
}

/*
 *  Construct z part of the evolution matrix
 */
void FloEvolIsingRandomSimpShiftCosRealTau::Evol_Z_Construct_(MatrixXcd & evol_z) {
    for (int i=0; i<dim_;i++){

        int state = i;
        int spin = 0, prev_spin = 0;
        double value = 0; // value which will be exponentiated

        for (int j=0; j<size_;j++){
            if (j>0) prev_spin = spin;
            spin = 2*(state & 1) - 1;
            state = state >> 1;

            value += random_h_[j] * spin;

            if (j>0) value += spin * prev_spin;
        }

        evol_z(i,i) = exp( - Complex_I * complex<double>(tau_ * value,0) );
    }
}








