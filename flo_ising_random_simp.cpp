//
// Created by Liangsheng Zhang on 4/14/15.
//

/*
 * Implementation of FloEvolIsingRandomSimp class
 */

#include <complex>
#include <iostream>
#include "constants.h"
#include "flo_evol_model.h"
#include "screen_output.h"
#include "randomc.h"

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

using namespace std;


/*
 * Construct the representation string and abstract type of the class.
 */
void FloEvolIsingRandomSimp::Repr_Init_(){
    repr_ << "Ising_Random_Simp_Floquet_L=" << size_ << ",W=" << W_;
    type_ = "Ising_Random_Simp_Floquet";
}

/*
 * Initialize random numbers
 */
void FloEvolIsingRandomSimp::Evol_Para_Init() {
    random_h_.resize(size_);

    // Random Fields
    for (int i=0; i<size_;i++){
        double u = RanGen_mersenne.Random();
        random_h_[i] = 1 + W_ * (2*u-1);
    }

    if (debug_){
        cout << "Random longitude field:" << endl;
        for (int i=0; i<size_;i++) cout << random_h_[i] << endl;
    }
}

/*
 * Copy random numbers
 */
void FloEvolIsingRandomSimp::Evol_Para_Copy(const vector< vector<double> >& random_h) {
    random_h_.resize(size_);

    cout << "Evol_Para_Copy not tested!!" << endl;
    abort();

    if(random_h.size() != 1){
        cout << type_ << " has only one set of random fields." << endl;
        cout << "Obtained: " << random_h.size() << endl;
        abort();
    }

    if(random_h[0].size() != size_){
        cout << "Incorrect number of random fields passed in." << endl;
        cout << "Expected: " << size_ << " Obtained: " << random_h[0].size() << endl;
        abort();
    }

    // Copy Random Fields and check the bounds
    const double lower_bound = 1 - W_;
    const double upper_bound = 1 + W_;
    for(int i=0; i<size_; i++){
        double val = random_h[0][i];
        if(val < lower_bound || val > upper_bound){
            cout << "Random field passed in at pos " << i << " is outside of bound." << endl;
            cout << "Lower bound: " << lower_bound << " Upper bound: " << upper_bound << endl;
            cout << "Val: " << val << endl;
            abort();
        }
        random_h_[i] = val;
    }

    if (debug_){
        cout << "Random longitude field:" << endl;
        for (int i=0; i<size_;i++) cout << random_h_[i] << endl;
    }
}

/*
 * Construct the time evolution matrix
 */
void FloEvolIsingRandomSimp::Evol_Construct() {

    if (!constructed_){
        evol_op_ = MatrixXcd::Zero(dim_, dim_);
        constructed_ = true;
    }
    else{
        cout << "Evolution operator has been constructed." << endl;
        abort();
    }

    MatrixXcd evol_x = MatrixXcd::Zero(dim_,dim_);
    MatrixXcd evol_z = MatrixXcd::Zero(dim_,dim_);

    Evol_X_Construct_(evol_x);
    Evol_Z_Construct_(evol_z);

    evol_op_ = evol_x * evol_z;

    if (debug_){
        cout << "X part of evolution matrix:" << endl;
        matrix_write(evol_x);
        cout << endl;

        cout << "Z part of evolution matrix:" << endl;
        matrix_write(evol_z);
        cout << endl;

        cout << "Evolution matrix:" << endl;
        matrix_write(evol_op_);
        cout << endl;
    }
}

/*
 * Construct X part of the evolution matrix
 */
void FloEvolIsingRandomSimp::Evol_X_Construct_(MatrixXcd & evol_x) {
    for (int i=0; i<dim_; i++){
        int flip_spin = 1;
        for (int j=0;j<size_;j++){
            int pos = i ^ flip_spin;
            evol_x(pos,i) += (1 - W_);
            flip_spin = flip_spin << 1;
        }
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<evol_x.cols(); i++){
        for (int j=i+1; j< evol_x.rows();j++){
            if (abs(evol_x(j,i) - evol_x(i,j)) > 1.0e-10){
                cout << "x part is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << evol_x(i,j) << endl;
                cout << "(j,i): " << evol_x(j,i) << endl;
                abort();
            }
        }
    }

    SelfAdjointEigenSolver<MatrixXcd> x_eigen;

    x_eigen.compute(evol_x);

    for (int i=0; i<dim_; i++){
        for (int j=0; j< dim_; j++){
            if (i==j){
                evol_x(i,j) = exp( - Complex_I * x_eigen.eigenvalues()[i] );
            }
            else{
                evol_x(i,j) = complex<double>(0,0);
            }
        }
    }

    evol_x = x_eigen.eigenvectors() * evol_x * x_eigen.eigenvectors().adjoint();
}

/*
 *  Construct z part of the evolution matrix
 */
void FloEvolIsingRandomSimp::Evol_Z_Construct_(MatrixXcd & evol_z) {
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

        evol_z(i,i) = exp( - Complex_I * complex<double>(value,0) );
    }
}

