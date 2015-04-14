//
// Created by Liangsheng Zhang on 4/14/15.
//


/*
 * This file implements some simple functions in flo_evol.h
 */

#include <iostream>
#include "flo_evol.h"

using namespace std;

void FloEvolVanilla::Transition_Compute(TransitionMatrix& transition,
                                        const string& matrix_name) const{
    if (matrix_name == "Basic_Full"){
        if (eigen_info_) transition.Basic_Full(eigen_ -> eigenvectors());
        else{
            cout << "Evolution Operator " << Repr() <<" has not been diagonalized with eigenvectors." << endl;
            abort();
        }
    }
    else{
        cout << "Transition matrix " << matrix_name << " cannot be established for "
             << Type() << endl;
        abort();
    }
}

void FloEvolVanilla::Evec(vector<MatrixXcd> evec) const {
    if (eigen_ == NULL){
        cout << Repr() << " has not been diagonalized with eigenvectors." << endl;
    }
    evec.resize(1);
    evec[0] = eigen_ -> eigenvectors();
}

void FloEvolVanilla::Eval(vector<VectorXcd> eval) const {
    if (eigen_ == NULL){
        cout << Repr() << " has not been diagonalized." << endl;
    }
    eval.resize(1);
    eval[0] = eigen_ -> eigenvalues();
}
