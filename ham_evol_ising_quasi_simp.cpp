//
// Created by Liangsheng Zhang on 9/29/15.
//

//
// Created by Liangsheng Zhang on 9/29/15.
//

/**
 ** Implementation of HamEvolIsingRandomSimpCos class
 **/

#include <complex>
#include <iostream>
#include "constants.h"
#include "ham_evol_model.h"
#include "screen_output.h"
#include "randomc.h"

using namespace std;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

/*
 * Construct the representation string and abstract type of the class.
 */
void HamEvolIsingQuasiSimp::Repr_Init_(){
    repr_ << "Ising_Quasi_Simp_Hamiltonian_L=" << size_ << ",W=" << W_;
    type_ = "Ising_Quasi_Simp_Hamiltonian";
}

/*
 * Initialize random fields
 */
void HamEvolIsingQuasiSimp::Evol_Para_Init() {
    random_h_.resize(size_);

    // Random phase
    double u = RanGen_mersenne.Random();
    double theta = 2 * Pi * u;

    // Random Fields
    for (int i=0; i<size_;i++){
        random_h_[i] = 1 + W_ * cos(i * Pi / Phi + theta);
    }

    if (debug_){
        cout << "Random phase: " << theta << endl;
        cout << "Longitude fields:" << endl;
        for (int i=0; i<size_; i++) cout << random_h_[i] << endl;
        cout << endl;
    }
}

/*
 * Construct the Hamiltonian
 */
void HamEvolIsingQuasiSimp::Evol_Construct() {
    if(!constructed_){
        ham_op_ = MatrixXd::Zero(dim_, dim_);
        constructed_ = true;
    }
    else{
        cout << "Evolution operator has been constructed." << endl;
        abort();
    }

    for(int i=0; i<dim_; i++){
        int state = i;
        int prev_spin = 0;
        int spin = 0;
        int x_state = 1; // Used to flip spin

        double z_value = 0; // diagonal element

        for (int j=0; j<size_;j++){
            if (j>0) prev_spin = spin;
            spin = 2*(state & 1) - 1;
            state = state >> 1;

            z_value += random_h_[j] * spin;
            if (j>0) z_value += spin * prev_spin;

            int flip_state = i ^ x_state; // Flip the spin
            x_state = x_state << 1;

            ham_op_(flip_state, i) += 1 - W_;
        }

        ham_op_(i, i) = z_value;
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<ham_op_.cols(); i++){
        for (int j=i+1; j< ham_op_.rows();j++){
            if (abs(ham_op_(j,i) - ham_op_(i,j)) > 1.0e-10){
                cout << "Hamiltonian is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << ham_op_(i,j) << endl;
                cout << "(j,i): " << ham_op_(j,i) << endl;
                abort();
            }
        }
    }

    if (debug_){
        cout << "Hamiltonian matrix:" << endl;
        matrix_write(ham_op_);
        cout << endl;
    }
}



