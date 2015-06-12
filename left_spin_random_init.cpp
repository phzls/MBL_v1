//
// Created by Liangsheng Zhang on 6/12/15.
//

#include <complex>
#include <iostream>
#include "init_obj.h"
#include "screen_output.h"

using namespace std;
using namespace Eigen;

/*
 * This function creates an initial state which has spin up on the leftmost spin and the
 * rest of the system is described by a density matrix proportional to identity matrix.
 */

void left_spin_random(const InitInfo& init_info, const TransitionMatrix& transition,
                      MatrixXcd& init_state_density, InitEvolData& init_evol_data){
    const int size = init_info.size; // System size
    const int total_rank = (1<<size); // Total dimension of Hilbert space
    const int rest_rank = (1<<(size-1)); // Total dimension of the system except the leftmost spin

    complex<double> diag = complex<double>(1.0/sqrt(rest_rank),0); // Diagonal elements

    // The initial density matrix written in basic binary basis
    MatrixXcd init_basic = MatrixXcd::Zero(total_rank, total_rank);

    for (int i=rest_rank; i<total_rank; i++) init_basic(i,i) = diag;

    // Change initial density matrix to evolution eigenvector basis
    init_state_density = transition.Matrix("Basic_Full").adjoint() * init_basic * transition.Matrix("Basic_Full");

    if ( (init_state_density.rows() != total_rank) || (init_state_density.cols() != total_rank)){
        cout << "Initial density matrix for left_spin_random has incorrect shape." << endl;
        cout << "Expected row/col num: " << total_rank << endl;
        cout << "Row num: " << init_state_density.rows() << endl;
        cout << "Col num: " << init_state_density.cols() << endl;
        abort();
    }

    for (int i=0; i< init_state_density.rows(); i++){
        for (int j=i; j<init_state_density.rows();j++){
            if ( norm( init_state_density(i,j) - conj(init_state_density(j,i)) ) > 1.0e-7 ){
                cout << "Initial density matrix for left_spin_random is not Hermitian at ("
                     << i << "," << j << ")." << endl;
                cout << "At (" << i << "," << j << "): " << init_state_density(i,j) << endl;
                cout << "At (" << j << "," << i << "): " <<  init_state_density(j,i) << endl;
                abort();
            }
        }
    }

    // Infinite time density matrix
    MatrixXcd diag_state = MatrixXcd::Zero(total_rank, total_rank);
    for (int i=0; i<total_rank; i++) diag_state(i,i) = init_state_density(i,i);

    // Reuse inital_basic
    init_basic = transition.Matrix("Basic_Full") * diag_state * transition.Matrix("Basic_Full").adjoint();

    init_evol_data.infinite_time_leftmost_spin = 0;

    for (int i=0; i<rest_rank;i++) {
        if (abs(imag(init_basic(i,i))) > 1.0e-7){
            cout << "Diagonal ensemble at " << i << "th diagonal in binary basis has imaginary part" << endl;
            cout << "Imaginary part: " << imag(init_basic(i,i)) << endl;
            abort();
        }
        init_evol_data.infinite_time_leftmost_spin -= real(init_basic(i,i));
    }
    for (int i=rest_rank; i<total_rank; i++) {
        if (abs(imag(init_basic(i,i))) > 1.0e-7){
            cout << "Diagonal ensemble at " << i << "th diagonal in binary basis has imaginary part" << endl;
            cout << "Imaginary part: " << imag(init_basic(i,i)) << endl;
            abort();
        }
        init_evol_data.infinite_time_leftmost_spin += real(init_basic(i,i));
    }

    if (init_info.debug){
        cout << "Initial state in density matrix:" << endl;
        complex_matrix_write(init_state_density);
        cout << endl;
        cout << "Infinite time leftmost spin:" << endl;
        cout << init_evol_data.infinite_time_leftmost_spin << endl;
        cout << endl;
    }

}