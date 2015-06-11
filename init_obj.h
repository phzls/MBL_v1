//
// Created by Liangsheng Zhang on 6/1/15.
//

#ifndef MBL_V1_INIT_OBJ_H
#define MBL_V1_INIT_OBJ_H

#include <map>
#include <Eigen/Dense>
#include "basis_transition.h"
#include "parameters.h"

using namespace std;
using namespace Eigen;

/**
 ** This file includes classes and functions used for constructing initial states.
 **/

/*
 * Information that is related to constructing initial states.
 */

struct InitInfo
{
    int size; // System size
    int dim; // Total dimension of Hilbert space
    double norm_delta; // A small number used to check whether norm is 1
    bool debug; // Whether output debug information
    string init_func_name; // The function name for initial state/density matrix construction

    // A deep copy constructor except for multi_ini_para
    InitInfo(const InitInfo&);

    // Default constructor
    InitInfo(): init_func_name("") {};

    // A deep copy operator
    void Copy(const InitInfo&);
};

// Pointer to all possible initial state construction function which gives a state vector
// in the basis of eigenstates
typedef void (*init_func)(const InitInfo&, const TransitionMatrix&, VectorXcd&);

// Pointer to all possible initial state construction function which gives a density matrix
// in the basis of eigenstates
typedef void (*init_func_C)(const InitInfo&, const TransitionMatrix&, MatrixXcd&);

class InitObj
{
private:
    // Map that stores init_func and its name to be called
    map<string, init_func>init_func_map_;

    // Map that stores init_func_C and its name to be called
    map<string, init_func_C>init_func_C_map_;

    void map_init_(); // Initialize init_func_map and init_func_C_map

public:
    InitInfo init_info;

    InitObj() {map_init_();}

    InitObj(const InitInfo& info): init_info(info) {map_init_();}

    // Construct initial state vector
    void Init_Func(const TransitionMatrix&, VectorXcd&) const;

    // Construct initial density matrix
    void Init_Func_C(const TransitionMatrix&, MatrixXcd&) const;

    // Print out all init_func
    void Print() const;

    // Print out all init_func_C
    void Print_C() const;

    ~InitObj(){};
};

/*
 * Initial state functions (more in old codes)
 */

void product_random(const InitInfo&, const TransitionMatrix&, VectorXcd&);

void random_product(const InitInfo&, const TransitionMatrix&, VectorXcd&);
void random_product(const InitInfo&, const TransitionMatrix&, MatrixXcd&);

void random_pure(const InitInfo&, const TransitionMatrix&, MatrixXcd&);

/*
 * Some useful functions for initial state construction
 */
// Construct a vector of amplitudes which can be used for random pure state
void random_pure_amplitude(vector<complex<double> >&);
// Check whehter the norm of a complex vector is close to 1
void norm_check(const VectorXcd&, double, const string&);
// Compute the complex density matrix of a corresponding normalized complex state vector
void state_to_density(const VectorXcd&, MatrixXcd&);


#endif //MBL_V1_INIT_OBJ_H

