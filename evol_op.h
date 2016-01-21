//
// Created by Liangsheng Zhang on 4/14/15.
//

#ifndef MBL_V1_EVOL_OP_H
#define MBL_V1_EVOL_OP_H

/*
 * A base class defines the evolution of a quantum system. The size of the system needs to
 * be passed in. In this case the local dimension at each site is assumed to be 2, and the
 * total dimension is calculated. The local dimension can also be specified when
 * constructing the class.
 */

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "basis_transition.h"

using namespace std;
using namespace Eigen;

class EvolOP
{
protected:
    const int size_; // Size of the system
    int dim_; // Dimension of the space given size
    const int local_dim_; // Local dimension of the Hilber space
    const bool iso_keep_; // Whether keep another EvolOP pointer which points
                          // to the isolated part of a model coupled to the bath

    // Pointers pointing to isolated part of a model coupling to the bath
    EvolOP* model_iso_;

    int model_id_;

public:
    // When local dimension is not given
    EvolOP(int size, bool iso_keep = false):
            size_(size), local_dim_(2), dim_(1 << size), iso_keep_(iso_keep),
            model_iso_(NULL) {model_num ++; model_id_ = model_num;}

    // When local dimension is explicitly given
    EvolOP(int size, int local_dim, bool iso_keep = false):
            size_(size), local_dim_(local_dim), dim_(int(pow(double(local_dim), double(size)))),
            iso_keep_(iso_keep), model_iso_(NULL) {model_num++; model_id_ = model_num;}

    vector<string> eigen_name; // Name of each eigensectors

    // Initialize parameters which will be used in constructing time evolution operator.
    // These parameters will only exist for concrete model classes.
    virtual void Evol_Para_Init() = 0;

    // Pass in a vector to initialize the matrix. It does not block the calling of Evol_Para_Init
    // which will then override the parameters
    virtual void Evol_Para_Copy(const vector< vector<double> >&) = 0;
    virtual void Evol_Para_Copy(const vector< vector< complex<double> > >&) = 0;

    // Constructing time evolution matrix
    virtual void Evol_Construct() = 0;

    // Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
    // False is not kept; True is kept. By default it is true.
    virtual void Evol_Diag(bool keep = true) = 0;

    // Compute various transition matrix. The string specifies what transition to compute
    virtual void Transition_Compute(TransitionMatrix&, const string&) const = 0;

    // Compute various transition matrix. The string specifies what transition to compute.
    // Eigenvectors are passed in, in case the internal eigenstates have been deleted.
    virtual void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXd>&) const = 0;

    // Compute various transition matrix. The string specifies what transition to compute
    // Eigenvectors are passed in, in case the internal eigenstates have been deleted.
    virtual void Transition_Compute(TransitionMatrix&, const string&, const vector<MatrixXcd>&) const = 0;

    // Erase the matrix to free some memroy
    virtual void OP_Erase() = 0;

    // Erase eigenvectors and eigenvalues
    virtual void Eigen_Erase() = 0;

    // Return the string format of representation string stream
    virtual string Repr() const = 0;

    // Return the type of the model as a string, i.e., representation without
    // concerete parameters
    virtual string Type() const = 0;

    // Return the type of eigenvectors: real or complex
    virtual string Eval_Type() const = 0;

    // Return the type of eigenvalues: real of complex
    virtual string Evec_Type() const = 0;

    // Return the type of basis which eigenstates are written in
    virtual string Eigen_Basis_Type() const = 0;

    // Return the size of the system
    int Get_Size() const {return size_;}

    // Return the dimension of total Hilbert space
    int Get_Dim() const {return dim_;}

    // Return the dimension of each sector of symmetry
    virtual vector<int> Get_Sector_Dim() const = 0;

    // Return eigenvectors
    virtual void Evec(vector<MatrixXcd>&) const = 0;
    virtual void Evec(vector<MatrixXd>&) const = 0;

    // Return eigenvalues
    virtual void Eval(vector<VectorXcd>&) const = 0;
    virtual void Eval(vector<VectorXd>&) const = 0;

    // Return the complex unitary operator according to a string
    virtual const MatrixXcd& Get_Complex(string) const = 0;

    // Return the complex unitary operator according to an integer
    virtual const MatrixXcd& Get_Complex(int) const = 0;

    // Return the real Hamiltonian operator according to a string
    virtual const MatrixXd& Get_Real(string) const = 0;

    // Return the real Hamiltonian operator according to an integer
    virtual const MatrixXd& Get_Real(int) const = 0;

    // Construct and return a Hamiltonian. The second string specifies the basis
    // and the last string specifies any extra requirement
    virtual void Get_Ham(MatrixXcd&, string, string) const {
        cout << "Get_Ham is not implemented for " << Repr() << endl;
        abort();
    }

    // Get the pointer to the complex isolated matrix
    void Get_Iso(EvolOP*& model) {model = model_iso_;}

    // Check whether eigensystems have been solved
    virtual bool Get_Eigen_Computed() const = 0;

    // Get the number for current model
    int Get_Model_ID() const {return model_id_;}

    static int model_num; // Model number

    virtual ~EvolOP(){
        if (model_iso_ != NULL) delete model_iso_;
    };
};


#endif //MBL_V1_EVOL_OP_H
