//
// Created by Liangsheng Zhang on 4/14/15.
//

#ifndef MBL_V1_FLO_EVOL_H
#define MBL_V1_FLO_EVOL_H

/*
 * This file contains base classes for Floquet time evolution operators
 */

#include <iostream>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "evol_op.h"

using namespace std;
using namespace Eigen;

/*
 * The evolution operator for Floquet system with no apparent symmetry to reduce
 * time evolution operator. It uses ComplexEigenSolver<MatrixXcd> for the public
 * member eigen inherited from EvolMatrix. It only uses eigen.
 */

class FloEvolVanilla : public EvolOP
{
protected:
    MatrixXcd evol_op_; // Time evolution operator
    bool constructed_; // Whether the matrix has been constructed and not erased

    stringstream repr_; // Representation string stream of the model
    string type_; // Type string of the model

    bool eigen_info_; // Whether eigenvectors have been computed

    ComplexEigenSolver<MatrixXcd>* eigen_;

public:
    // When local dimension is not given
    FloEvolVanilla(int size): EvolOP(size), constructed_(false), eigen_info_(false), eigen_(NULL)
    {eigen_name.resize(1,"");}

    // When local dimension is given
    FloEvolVanilla(int size, int local_dim): EvolOP(size, local_dim), constructed_(false),
                                             eigen_info_(false), eigen_(NULL) {eigen_name.resize(1,"");}

    // Diagnolize time evolution matrix with eigenvectors kept
    void Evol_Diag(bool keep = true) {
        if (constructed_){
            if (eigen_ == NULL){
                eigen_ = new ComplexEigenSolver<MatrixXcd>;
                eigen_ -> compute(evol_op_, keep);
                eigen_info_ = keep;
                eigen_name[0] = "Full";
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

    // Return the string format of representation string stream
    string Repr() const {return repr_.str();}

    // Return the type of the model
    string Type() const {return type_;}

    // Return the type of the basis that eigenstates are written in
    string Eigen_Type() const {return "Basic";}

    // Return dimension of each sector. Here only 1 sector exists, so total dimension
    // is returned.
    vector<int> Get_Sector_Dim() const{
        vector<int> dim(1);
        dim[0] = dim_;
        return dim;
    }

    // Erase the evolutionary operator
    void OP_Erase() {evol_op_.resize(0,0); constructed_ = false;}

    // Erase eigenvectors and eigenvalues
    void Eigen_Erase() { if (eigen_ != NULL) {delete eigen_; eigen_ = NULL;}}

    // Construct Transition Matrix
    void Transition_Compute(TransitionMatrix&, const string&) const;

    // Return eigenvectors in a column_wise fashion at each index
    void Evec(vector<MatrixXcd>) const;
    void Evec(vector<MatrixXd> evec) const{
        cout << Repr() << " has no real eigenvectors." << endl;
        abort();
    }

    // Return eigenvalues
    void Eval(vector<VectorXcd>) const;
    void Eval(vector<VectorXd>) const{
        cout << Repr() << " has no real eigenvalues." << endl;
        abort();
    }

    // Return evol_op no matter what the input is
    const MatrixXcd& Get_Complex(string a) const {
        if (constructed_) return evol_op_;
        else{
            cout << Repr() << " has not been constructed yet." << endl;
            abort();
        }
    }

    const MatrixXcd& Get_Complex(int a=0) const {
        if (constructed_) return evol_op_;
        else{
            cout << Repr() << " has not been constructed yet." << endl;
            abort();
        }
    }

    // No real Hamiltonian to return
    const MatrixXd& Get_Real(string a) const{
        cout << "No real operator is constructed for " << Repr() << endl;
        abort();
    }

    const MatrixXd& Get_Real(int a = 0) const{
        cout << "No real operator is constructed for " << Repr() << endl;
        abort();
    }

    virtual ~FloEvolVanilla() { if (eigen_ != NULL) {delete eigen_; eigen_ = NULL;}}
};

#endif //MBL_V1_FLO_EVOL_H
