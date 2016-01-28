//
// Created by Liangsheng Zhang on 4/15/15.
//

#ifndef MBL_V1_METHODS_H
#define MBL_V1_METHODS_H

/*
 * This file contains methods which are used in various main tasks.
 */

#include <utility>
#include <vector>
#include <Eigen/Core>
#include "evol_op.h"
#include "parameters.h"

using namespace std;
using namespace Eigen;


/*
 * Write eigenstates of an evolutionary operator in basic binary basis, and output the result
 * in passed in vectors. The components of each eigenvector is stored in the inner index
 */
void evec_to_basic(const EvolOP*, const vector<MatrixXd>&, vector<vector<double> >&);
void evec_to_basic(const EvolOP*, const vector<MatrixXd>&, vector<vector<complex<double> > >&);
void evec_to_basic(const EvolOP*, const vector<MatrixXcd>&, vector<vector<complex<double> > >&);

/*
 * Construct a left reduced density matrix from a state vector for a spin chain where local
 * dimension is 2. The first integer is the total size of spin chain, and the second integer
 * is the size of left part.
 */
void reduced_density_left_2(const VectorXcd&, int, int, MatrixXcd&);
void reduced_density_left_2(const VectorXd&, int, int, MatrixXcd&);
void reduced_density_left_2(const VectorXd&, int, int, MatrixXd&);

/*
 * Construct a vector of left-side binary entropies for a given model of spin chain where
 * local dimension is 2. The integer is the size of the left part.
 */
void model_entropy_left_2(const EvolOP*, int, vector<double>&);

/*
 * Construct the left entropy of a given vector state from a model of spin chain where
 * local dimension is 2. The first integer is the total size of the spin chain, and
 * the second is the size of left part.
 */
double state_entropy_left_2(const VectorXd&, int, int);
double state_entropy_left_2(const VectorXcd&,int, int);

/*
 * Generating random parameters used for various models according to the model name
 * and possible configuration parameters. The number of parameters generated for
 * each inner vector is given by the integer
 */
void model_para_generation(const AllPara&, vector< vector<double> >&, int);


#endif //MBL_V1_METHODS_H
