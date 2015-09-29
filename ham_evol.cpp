//
// Created by Liangsheng Zhang on 9/17/15.
//

/*
 * This file implements some simple functions in flo_evol.h
 */

#include <iostream>
#include "ham_evol.h"

using namespace std;

void HamEvolVanillaReal::Evec(vector<MatrixXd>& evec) const {
    if (eigen_ == NULL){
        cout << Repr() << " has not been diagonalized with eigenvectors." << endl;
        abort();
    }
    evec.resize(1);
    evec[0] = eigen_ -> eigenvectors();
}

void HamEvolVanillaReal::Eval(vector<VectorXd>& eval) const {
    if (eigen_ == NULL){
        cout << Repr() << " has not been diagonalized." << endl;
        abort();
    }
    eval.resize(1);
    eval[0] = eigen_ -> eigenvalues();
}

void HamEvolVanillaReal::Transition_Compute(TransitionMatrix& transition, const string& matrix_name) const{
    if (matrix_name == "Basic_Full"){
        if (eigen_info_) transition.Basic_Full(eigen_ -> eigenvectors());
        else{
            cout << "Evolution Operator " << Repr() <<" has not been diagonalized." << endl;
            abort();
        }
    }
    else{
        cout << "Transition matrix " << matrix_name << " cannot be established for "
        << Type() << endl;
        abort();
    }
}

void HamEvolVanillaReal::Transition_Compute(TransitionMatrix& transition, const string& matrix_name,
                                            const vector<MatrixXcd>& evec) const{
    if (matrix_name == "Basic_Full"){
        if (eigen_info_) transition.Basic_Full(eigen_ -> eigenvectors());
        else transition.Basic_Full(evec[0]);
    }
    else{
        cout << "Transition matrix " << matrix_name << " cannot be established for "
        << Type() << endl;
        abort();
    }
}

void HamEvolVanillaReal::Transition_Compute(TransitionMatrix& transition, const string& matrix_name,
                                            const vector<MatrixXd>& evec) const{
    if (matrix_name == "Basic_Full"){
        if (eigen_info_) transition.Basic_Full(eigen_ -> eigenvectors());
        else transition.Basic_Full(evec[0]);
    }
    else{
        cout << "Transition matrix " << matrix_name << " cannot be established for "
        << Type() << endl;
        abort();
    }
}
