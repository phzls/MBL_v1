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

void FloEvolVanilla::Transition_Compute(TransitionMatrix& transition, const string& matrix_name,
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

void FloEvolVanilla::Transition_Compute(TransitionMatrix& transition, const string& matrix_name,
                                        const vector<MatrixXd>& evec) const{
    cout << "Eigenvectors for transition computation of " << Type() << " must be complex." << endl;
    abort();
}

void FloEvolVanilla::Evec(vector<MatrixXcd>& evec) const {
    if (eigen_ == NULL){
        cout << Repr() << " has not been diagonalized with eigenvectors." << endl;
    }
    evec.resize(1);
    evec[0] = eigen_ -> eigenvectors();
}

void FloEvolVanilla::Eval(vector<VectorXcd>& eval) const {
    if (eigen_ == NULL){
        cout << Repr() << " has not been diagonalized." << endl;
    }
    eval.resize(1);
    eval[0] = eigen_ -> eigenvalues();
}




//=========================================================================================





void FloEvolVanillaReal::Transition_Compute(TransitionMatrix& transition, const string& matrix_name) const{
    if (matrix_name == "Basic_Full"){
        if (eigen_info_) transition.Basic_Full(evec_);
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

void FloEvolVanillaReal::Transition_Compute(TransitionMatrix& transition, const string& matrix_name,
                                            const vector<MatrixXcd>& evec) const{
    if (matrix_name == "Basic_Full"){
        if (eigen_info_) transition.Basic_Full(evec_);
        else transition.Basic_Full(evec[0]);
    }
    else{
        cout << "Transition matrix " << matrix_name << " cannot be established for "
        << Type() << endl;
        abort();
    }
}

void FloEvolVanillaReal::Transition_Compute(TransitionMatrix& transition, const string& matrix_name,
                                            const vector<MatrixXd>& evec) const{
    if (matrix_name == "Basic_Full"){
        if (eigen_info_) transition.Basic_Full(evec_);
        else transition.Basic_Full(evec[0]);
    }
    else{
        cout << "Transition matrix " << matrix_name << " cannot be established for "
        << Type() << endl;
        abort();
    }
}

void FloEvolVanillaReal::Evol_Diag(bool keep){
    if (constructed_){
        if (!diag_){
            SelfAdjointEigenSolver<MatrixXd>* eigen_; // Only one sector

            eigen_ = new SelfAdjointEigenSolver<MatrixXd>;
            eigen_ -> compute(evol_op_real_);

            eigen_info_ = keep;
            eigen_name[0] = "Full";

            eval_.resize(dim_);
            VectorXd temp(dim_);

            for (int i=0; i < eval_.rows(); i++){
                double real = eigen_ -> eigenvalues()[i];

                temp = evol_op_imag_ * eigen_ -> eigenvectors().col(i);

                int index = 0;
                while (abs(temp[index]) < 2.0e-6 ) index ++;

                double imag = temp[index] / ( eigen_ -> eigenvectors().col(i)[index] );

                eval_[i] = complex<double>(real, imag);
            }

            if (keep){
                evec_.resize(dim_, dim_);
                for (int i=0; i<dim_;i++){
                    for (int j=0; j<dim_; j++){
                        evec_(j,i) = eigen_ -> eigenvectors()(j,i);
                    }
                }
            }

            delete eigen_;
            eigen_ = NULL;
            diag_ = true;
        }
        else{
            cout << "Matrix has already been diagonalized." << endl;
            abort();
        }
    }
    else{
        cout << "The matrix for diagonalization does not exist." <<endl;
        abort();
    }
}

const MatrixXd& FloEvolVanillaReal::Get_Real(int num) const {
    if (constructed_){
        if (num == 0){
            return evol_op_real_;
        }
        else if (num == 1){
            return evol_op_imag_;
        }
        else{
            cout << "The type of time evolution operator is not understood for "
            << Repr() << endl;
            abort();
        }
    }
    else{
        cout << Repr() << " has not been constructed yet." << endl;
        abort();
    }
}

const MatrixXd& FloEvolVanillaReal::Get_Real(string name) const {
    if (constructed_){
        if (name == "real" || name == "Real"){
            return evol_op_real_;
        }
        else if (name == "imag" || name == "Imag"){
            return evol_op_imag_;
        }
        else{
            cout << "The type of time evolution operator is not understood." << endl;
            abort();
        }
    }
    else{
        cout << Repr() << " has not been constructed yet." << endl;
        abort();
    }
}

void FloEvolVanillaReal::Evec(vector<MatrixXd>& evec) const {
    if (!eigen_info_){
        cout << Repr() << " has not been diagonalized with eigenvectors." << endl;
    }
    evec.resize(1);
    evec[0] = evec_;
}

void FloEvolVanillaReal::Evec(vector<MatrixXcd>& evec) const {
    if (!eigen_info_){
        cout << Repr() << " has not been diagonalized with eigenvectors." << endl;
    }
    evec.resize(1);
    evec[0].resize(evec_.rows(), evec_.cols());

    for (int i=0; i<evec_.cols(); i++){
        for (int j=0; j<evec_.rows(); j++){
            evec[0](j,i) = complex<double>(evec_(j,i),0);
        }
    }

}

void FloEvolVanillaReal::Eval(vector<VectorXcd>& eval) const {
    if (!diag_){
        cout << Repr() << " has not been diagonalized." << endl;
    }
    eval.resize(1);
    eval[0] = eval_;
}
