//
// Created by Liangsheng Zhang on 12/3/15.
//

//
// Created by Liangsheng Zhang on 12/1/15.
//

/*
 * Implementation of FloEvolXXZUniformRandomShiftReal class
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
void FloEvolXXZUniformRandomShiftReal::Repr_Init_(){
    repr_ << "XXZ_Uniform_Random_Shift_Real_Floquet_L=" << size_ << ",J=" << J_;
    type_ = "XXZ_Uniform_Random_Shift_Real_Floquet";
}

/*
 * Initialize random numbers
 */
void FloEvolXXZUniformRandomShiftReal::Evol_Para_Init() {
    uniform_.resize(size_);

    // Random Fields
    for (int i=0; i<size_; i++){
        double u = RanGen_mersenne.Random();
        uniform_[i] = 2*sqrt(3)*u - sqrt(3);
    }

    J_ = sqrt( 1 - sigma_*sigma_ * ( (h_*h_ + g_*g_)/( 1 - 1.0/size_ + g_*g_ + h_*h_ ) ) );

    if (debug_){
        cout << "Uniform random numbers:" << endl;
        for (int i=0; i<size_;i++) cout << uniform_[i] << endl;
        cout << "sigma: " << sigma_ << endl;
        cout << "J: " << J_ << endl;
    }

    initialized_ = true;
}

/*
 * Construct the time evolution matrix
 */
void FloEvolXXZUniformRandomShiftReal::Evol_Construct() {
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
void FloEvolXXZUniformRandomShiftReal::Evol_X_Construct_(MatrixXcd & evol_half_x) {
    const int int_size = sizeof(int)*CHAR_BIT; // Number of bits in an integer
    for(int i=0; i<dim_; i++){
        bitset<int_size> spin1(i);

        for(int j=0; j<dim_; j++){
            bitset<int_size> spin2(j);
            complex<double> prod = 1;

            for(int k=0; k<size_; k++){
                if(spin1[k] == spin2[k]) prod *= cos( tau_*g_*(J_+sigma_*uniform_[k])/2.0 );
                else prod *= complex<double>( 0, -sin( tau_*g_*(J_+sigma_*uniform_[k])/2.0 ) );
            }

            evol_half_x(i,j) = prod;
            if(i!=j) evol_half_x(j,i) = prod;
        }
    }
}

/*
 *  Construct z part of the evolution matrix
 */
void FloEvolXXZUniformRandomShiftReal::Evol_Z_Construct_(MatrixXcd & evol_z) {
    const int int_size = sizeof(int)*CHAR_BIT; // Number of bits in an integer
    for (int i=0; i<dim_;i++){
        bitset<int_size> state(i);
        double value = 0;
        for (int j=0; j<size_;j++){
            value += h_*(J_+sigma_*uniform_[j])*(2*state[j]-1);
            if(j>0) value += J_*(2*state[j]-1)*(2*state[j-1]-1);
        }

        evol_z(i,i) = exp( complex<double>(0,-tau_*value) );
    }
}

/*
 * Construct the following Hamiltonian no matter what the input string is
 */
void FloEvolXXZUniformRandomShiftReal::Get_Ham(MatrixXcd& ham, string basis, string s){
    ham = MatrixXcd::Zero(dim_, dim_);
    const int int_size = sizeof(int)*CHAR_BIT; // Number of bits in an integer

    // Construct the matrix in basic basis
    for(int i=0; i<dim_; i++){
        bitset<int_size> state(i);

        double z_value = 0; // diagonal element

        for (int j=0; j<size_;j++){
            z_value += h_*(J_+sigma_*uniform_[j])*(2*state[j]-1);
            if(j>0) z_value += J_*(2*state[j]-1)*(2*state[j-1]-1);

            state.flip(j);
            int x_state = state.to_ulong();
            state.flip(j); // Flip back the spin

            ham(x_state, i) += g_*(J_+sigma_*uniform_[j]);
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



