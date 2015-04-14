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

using namespace std;
using namespace Eigen;

class EvolOP
{
protected:
    const int size_; // Size of the system
    const int dim_; // Dimension of the space given size
    const int local_dim_; // Local dimension of the Hilber space
    const bool iso_keep_; // Whether keep another EvolOP pointer which points
                          // to the isolated part of a model coupled to the bath
    bool eigen_computed_; // Whether eigensystems of this model have been solved

    // Pointers pointing to isolated part of a model coupling to the bath
    EvolOP* model_iso_;

public:
    // When local dimension is not given
    EvolOP(int size, bool iso_keep = false):
            size_(size), local_dim_(2), dim_(1 << size), iso_keep_(iso_keep),
            model_iso_(NULL),eigen_computed_(false) {};

    // When local dimension is explicitly given
    EvolOP(int size, int local_dim, bool iso_keep = false):
            size_(size), local_dim_(local_dim), dim_(int(pow(double(local_dim), double(size)))),
            iso_keep_(iso_keep), model_iso_(NULL), eigen_computed_(false) {};

    vector<string> eigen_name; // Name of each eigensectors

    // Initialize parameters which will be used in constructing time evolution operator.
    // These parameters will only exist for concrete model classes.
    virtual void Evol_Para_Init() = 0;

    // Constructing time evolution matrix
    virtual void Evol_Construct() = 0;

    // Diagnolize time evolution matrix, user can determine whether eigenvectors are kept
    // False is not kept; True is kept. By default it is true.
    virtual void Evol_Diag(bool keep = true) = 0;

    // Compute various transition matrix. The string specifies what transition to compute
    virtual void Transition_Compute(TransitionMatrix&, const string&) const = 0;

    // Erase the matrix to free some memroy
    virtual void Evol_Erase() = 0;

    // Return the string format of representation string stream
    virtual string Repr() const = 0;

    // Return the type of the model as a string, i.e., representation without
    // concerete parameters
    virtual string Type() const = 0;

    // Return the type of basis which eigenstates are written in
    virtual string Eigen_Type() const = 0;

    // Return the size of the system
    int Get_Size() const {return size_;}

    // Return the dimension of total Hilbert space
    int Get_Dim() const {return dim_;}

    // Return the dimension of each sector of symmetry
    virtual vector<int> Get_Sector_Dim() const = 0;

    // Return the complex unitary operator according to a string
    virtual const MatrixXcd& Get_Complex(string) const = 0;

    // Return the complex unitary operator according to an integer
    virtual const MatrixXcd& Get_Complex(int) const = 0;

    // Return the real Hamiltonian operator according to a string
    virtual const MatrixXd& Get_Real(string) const = 0;

    // Return the real Hamiltonian operator according to an integer
    virtual const MatrixXd& Get_Real(int) const = 0;

    // Get the pointer to the complex isolated matrix
    void Get_Iso(EvolOP*& model) {model = model_iso_;}

    // Check whether eigensystems have been solved
    bool Get_Eigen_Computed() const {return eigen_computed_;}

    virtual ~EvolOP(){
        if (model_iso_ != NULL) delete model_iso_;
    };
};


#endif //MBL_V1_EVOL_OP_H
