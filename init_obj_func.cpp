//
// Created by Liangsheng Zhang on 6/1/15.
//

/**
 ** This file contains which are useful for initial state construction.
 **/

#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include "constants.h"
#include "randomc.h"
#include "init_obj.h"

using namespace std;
using namespace Eigen;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

/*
 * This function converts a complex state vector to a complex density matrix
 */

void state_to_density(const VectorXcd& state, MatrixXcd& matrix){
    matrix = MatrixXcd::Zero(state.size(), state.size());

    for (int i=0; i<state.size(); i++){
        for (int j=i; j< state.size(); j++){
            matrix(i,j) = state(i) * conj(state(j));
            if (i!=j) matrix(j,i) = conj(matrix(i,j));
        }
    }
}

/*
 * This function checks whether total norm of a complex vector is close 1. delta gives the small
 * number for individual component error tolerance. state_type gives the name of the state.
 */

void norm_check(const VectorXcd& state, double delta, const string& state_type){
    double state_norm = 0;
    for (int i=0; i< state.size(); i++){
        state_norm += norm(state(i));
    }

    cout << state_type << " norm: " << state_norm << endl;

    // Consider delta to be the error tolerance for the individual component
    if (abs(state_norm - 1) > (delta * state.size())){
        cout << state_type << " norm is not 1." << endl;
        abort();
    }
}

/*
 * This function creates a vector of random complex amplitudes which can be used for random
 * pure state
 */

void random_pure_amplitude(vector<complex<double> >& amplitude){
    double sum = 0;
    for (int i=0;i<amplitude.size();i++)
    {
        double U1 = RanGen_mersenne.Random();
        double U2 = RanGen_mersenne.Random();

        double real = sqrt(-2*log(1-U1))*cos(2*Pi*U2);
        double imag = sqrt(-2*log(1-U1))*sin(2*Pi*U2);

        amplitude[i] = complex<double>(real,imag);

        sum += norm(amplitude[i]);
    }

    for (int i=0;i<amplitude.size();i++)
    {
        amplitude[i] /= complex<double>(sqrt(sum),0);
    }
}

