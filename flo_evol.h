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

    ComplexEigenSolver<MatrixXcd>* eigen_; // Only one sector

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

    // Return the type of eigenvalues
    string Eval_Type() const {return "Complex";}

    // Return the type of eigenvectors
    string Evec_Type() const {return "Complex";}

    // Return the type of the basis that eigenstates are written in
    string Eigen_Basis_Type() const {return "Basic";}

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
    void Eigen_Erase() { if (eigen_ != NULL) {delete eigen_; eigen_ = NULL; eigen_info_ = false;}}

    // Construct Transition Matrix
    void Transition_Compute(TransitionMatrix&, const string&) const;
    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXcd>&) const;
    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXd>&) const;

    // Return eigenvectors in a column_wise fashion at each index
    void Evec(vector<MatrixXcd>&) const;
    void Evec(vector<MatrixXd>& evec) const{
        cout << Repr() << " has no real eigenvectors." << endl;
        abort();
    }

    // Return eigenvalues
    void Eval(vector<VectorXcd>&) const;
    void Eval(vector<VectorXd>& eval) const{
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
        cout << "No real Hamiltonian is constructed for " << Repr() << endl;
        abort();
    }

    const MatrixXd& Get_Real(int a = 0) const{
        cout << "No real Hamiltonian is constructed for " << Repr() << endl;
        abort();
    }

    bool Get_Eigen_Computed() const {return (eigen_ != NULL);}

    virtual ~FloEvolVanilla() { if (eigen_ != NULL) {delete eigen_; eigen_ = NULL;}}
};

#endif //MBL_V1_FLO_EVOL_H





//===================================================================================================





/*
 * The evolution operator for Floquet system with no apparent symmetry to reduce
 * time evolution operator. Its eigenvectors are real, so it can be diagnalized
 * using only the real part of the evolution matrix. The real and imaginary parts
 * of it must be Hermitian matrix. It uses SelfAdjointEigenSolver<MatrixXdd> for the public
 * member eigen inherited from EvolMatrix. It only uses eigen. eigenvalues in
 * eigen member will only have real part. The imaginary part can be computed
 * from the imaginary part of the operator. Both real and imaginary parts
 * can be accessed from Get_Real_H() function. All three matrices will be
 * preserved until erased. 0 is for real part of the matrix, 1 for imaginary part.
 */

class FloEvolVanillaReal : public EvolOP
{
protected:
    MatrixXcd evol_op_; // Time evolution operator
    MatrixXd evol_op_real_; // Real part of time evolution operator
    MatrixXd evol_op_imag_; // Imaginary part of time evolution operator

    VectorXcd eval_; // Eigenvalues
    MatrixXd evec_; // Eigenvectors

    bool constructed_; // Whether the matrix has been constructed and not erased
    bool diag_; // Whether it is diagonalized

    stringstream repr_; // Representation string stream of the model
    string type_; // Type string of the model

    bool eigen_info_; // Whether eigenvectors have been computed

public:
    // When local dimension is not given
    FloEvolVanillaReal(int size): EvolOP(size), constructed_(false), eigen_info_(false)
    {eigen_name.resize(1,"");}

    // When local dimension is given
    FloEvolVanillaReal(int size, int local_dim): EvolOP(size, local_dim), constructed_(false), eigen_info_(false)
    {eigen_name.resize(1,"");}

    // Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
    // False is not kept; True is kept
    void Evol_Diag(bool keep = true);

    // Return the string format of representation string stream
    string Repr() const {return repr_.str();}

    // Return the type of the model
    string Type() const {return type_;}

    // Return the type of eigenvalues
    string Eval_Type() const {return "Complex";}

    // Return the type of eigenvectors
    string Evec_Type() const {return "Real";}

    // Return the type of the basis that eigenstates are written in
    string Eigen_Basis_Type() const {return "Basic";}

    // Return dimension of each sector. Here only 1 sector exists, so total dimension
    // is returned.
    vector<int> Get_Sector_Dim() const{
        vector<int> dim(1);
        dim[0] = dim_;
        return dim;
    }

    // Erase the evolutionary operator and its real and imaginary parts
    void OP_Erase() {evol_op_.resize(0,0); evol_op_real_.resize(0,0);
        evol_op_imag_.resize(0,0); constructed_ = false;}

    // Erase eigenvectors and eigenvalues
    void Eigen_Erase() { eval_.resize(0); evec_.resize(0,0); diag_ = false; eigen_info_ = false;}

    // Construct Transition Matrix
    void Transition_Compute(TransitionMatrix&, const string&) const;
    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXd>&) const;
    void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXcd>&) const;

    // Return eigenvectors in a column_wise fashion at each index
    void Evec(vector<MatrixXd>&) const;
    void Evec(vector<MatrixXcd>&) const;

    // Return eigenvalues
    void Eval(vector<VectorXcd>&) const;
    void Eval(vector<VectorXd>& eval) const{
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

    bool Get_Eigen_Computed() const {return diag_;}

    const MatrixXd& Get_Real(string a) const;

    const MatrixXd& Get_Real(int a = 0) const;

    virtual ~FloEvolVanillaReal() {};
};
