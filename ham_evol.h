//
// Created by Liangsheng Zhang on 9/17/15.
//

#ifndef MBL_V1_HAM_EVOL_H
#define MBL_V1_HAM_EVOL_H

/*
 * This file contains base classes for Hamiltonian time evolution operators
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
 * The evolution operator for real Hamiltonian system with no apparent symmetry to reduce
 * time evolution operator, so it will have only one sector: 0 sector. THe Hamiltonian is real,
 * so its eigenvectors are real, It uses SelfAdjointEigenSolver<MatrixXd> for the public member
 * eigen inherited from EvolMatrix. It only uses eigen. eigenvalues in eigen member will only have
 * real part. Hamiltonian can be accessed from Get_Real_H() function. All matrices will be
 * preserved until erased.
 */

class HamEvolVanillaReal : public EvolOP
{
protected:
    MatrixXd ham_op_; // Hamiltonian operator

    SelfAdjointEigenSolver<MatrixXd>* eigen_; // Only one sector

    bool constructed_; // Whether the matrix has been constructed and not erased

    stringstream repr_; // Representation string stream of the model
    string type_; // Type string of the model

    bool eigen_info_; // Whether eigenvectors have been computed

public:
    // When local dimension is not given
    HamEvolVanillaReal(int size): EvolOP(size), constructed_(false), eigen_info_(false)
    {eigen_name.resize(1,"");}

    // When local dimension is given
    HamEvolVanillaReal(int size, int local_dim): EvolOP(size, local_dim), constructed_(false), eigen_info_(false)
                                                 {eigen_name.resize(1,"");}

    // Diagnolize Hamiltonian matrix, user can determine whether eigenvectors are kept
    // False is not kept; True is kept
    void Evol_Diag(bool keep = true) {
        if (constructed_){
            if (eigen_ == NULL){
                eigen_ = new SelfAdjointEigenSolver<MatrixXd>;
                int evec_flag = Eigen::EigenvaluesOnly;
                if(keep) evec_flag = Eigen::ComputeEigenvectors;

                eigen_ -> compute(ham_op_, evec_flag);
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

    // Return the type of eigenvalues
    string Eval_Type() const {return "Real";}

    // Return the type of eigenvectors
    string Evec_Type() const {return "Real";}

    // Return dimension of each sector. Here only 1 sector exists, so total dimension
    // is returned.
    vector<int> Get_Sector_Dim() const{
        vector<int> dim(1);
        dim[0] = dim_;
        return dim;
    }

    // Erase the evolutionary operator
    void OP_Erase() {ham_op_.resize(0,0); constructed_ = false;}

    // Erase eigenvectors and eigenvalues
    void Eigen_Erase() { if (eigen_ != NULL) {delete eigen_; eigen_ = NULL; eigen_info_ = false;}}

    // Return eigenvectors in a column_wise fashion at each index
    void Evec(vector<MatrixXd>&) const;
    void Evec(vector<MatrixXcd>&) const{
        cout << Repr() << "has no complex eigenvectors." << endl;
        abort();
    }

    // Return eigenvalues
    void Eval(vector<VectorXd>&) const;
    void Eval(vector<VectorXcd>& eval) const{
        cout << Repr() << " has no complex eigenvalues." << endl;
        abort();
    }

    // Return ham_op no matter what the input is
    const MatrixXd& Get_Real(string a) const {
        if (constructed_) return ham_op_;
        else{
            cout << Repr() << " has not been constructed yet." << endl;
            abort();
        }
    }

    const MatrixXd& Get_Real(int a=0) const {
        if (constructed_) return ham_op_;
        else{
            cout << Repr() << " has not been constructed yet." << endl;
            abort();
        }
    }

    bool Get_Eigen_Computed() const {return (eigen_ != NULL);}

    const MatrixXcd& Get_Complex(string a) const {
        cout << "No complex Hamiltonian is constructed for " << Repr() << endl;
        abort();
    }

    const MatrixXcd& Get_Complex(int a = 0) const {
        cout << "No complex Hamiltonian is constructed for " << Repr() << endl;
        abort();
    }

    virtual ~HamEvolVanillaReal() { if (eigen_ != NULL) {delete eigen_; eigen_ = NULL;} };
};


#endif //MBL_V1_HAM_EVOL_H
