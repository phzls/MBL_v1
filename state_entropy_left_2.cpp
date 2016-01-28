//
// Created by Liangsheng Zhang on 1/27/16.
//

/*
 * This file calculates the entropy for the left part of a spin chain from a given vector,
 * where the spin chain has local dimension at each site 2. Size is the total chain length
 * and left_size is the length of the left part. It assumes the passed in vector is in a
 * basic binary representation. This representation has basis states which have 0 being down
 * and 1 being up, and the counting starts from the rightmost spin.
 */

#include <iostream>
#include "methods.h"

using namespace std;
using namespace Eigen;

double state_entropy_left_2(const VectorXd& state, int size, int left_size){
    MatrixXd reduced_d; // Reduced density matrix
    reduced_density_left_2(state, size, left_size, reduced_d);

    double ent = 0; // Entropy

    SelfAdjointEigenSolver<MatrixXd> density_eigen; // Eigen for reduced density matrix
    density_eigen.compute(reduced_d, false); // Eigenvectors not computed

    // Compute entropy for this state
    for (int j=0; j<density_eigen.eigenvalues().rows();j++){
        double eval = density_eigen.eigenvalues()(j);
        if (abs(eval)>1.0e-10)
        {
            if (eval<0){
                cout << "Density matrix has significant negative eigenvalues." << endl;
                cout << eval << endl;
                abort();
            }
            ent += -eval*log2(eval);
        }
    }

    return ent;
}

double state_entropy_left_2(const VectorXcd& state, int size, int left_size){
    MatrixXcd reduced_d; // Reduced density matrix
    reduced_density_left_2(state, size, left_size, reduced_d);

    double ent = 0; // Entropy

    SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix
    density_eigen.compute(reduced_d, false); // Eigenvectors not computed

    // Compute entropy for this state
    for (int j=0; j<density_eigen.eigenvalues().rows();j++){
        double eval = density_eigen.eigenvalues()(j);
        if (abs(eval)>1.0e-10)
        {
            if (eval<0){
                cout << "Density matrix has significant negative eigenvalues." << endl;
                cout << eval << endl;
                abort();
            }
            ent += -eval*log2(eval);
        }
    }

    return ent;
}

