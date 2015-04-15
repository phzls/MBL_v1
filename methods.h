//
// Created by Liangsheng Zhang on 4/15/15.
//

#ifndef MBL_V1_METHODS_H
#define MBL_V1_METHODS_H

/*
 * This file contains methods which are used in evolution/level statistics.
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

#endif //MBL_V1_METHODS_H
