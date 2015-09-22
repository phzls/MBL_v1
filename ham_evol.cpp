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
