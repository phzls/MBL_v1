//
// Created by Liangsheng Zhang on 4/15/15.
//

#ifndef MBL_V1_GENERIC_FUNC_H
#define MBL_V1_GENERIC_FUNC_H

/*
 * This file contains some small functions which generalize some functions from c++ so that
 * they can take a larger range of arguments.
 */

#include <complex>
#include <cmath>
#include <vector>

using namespace std;

/*
 * Generic conjugate function which also handles double and integer
 */
template<class T>
complex<T> generic_conj(complex<T> val){ return conj(val);}

double generic_conj(double);
int generic_conj(int);

/*
 * Generic norm function which also handles double and integer. Note the norm in c++
 * is the square of the unsual definition of norm, so here abs is used.
 */
template<class T>
T generic_norm(complex<T> val){ return abs(val);}

double generic_norm(double);
int generic_norm(int);

/*
 * Given a vector of data, compute its sample mean and the standard deviation from sample variance
 */
template<class T>
void generic_mean_sd(const vector<T>& data, T& mean, T& sd){
    mean = 0;
    sd = 0;

    const int N = data.size();

    for (int i=0; i<N; i++){
        mean += data[i];
        sd += data[i] * data[i];
    }

    mean /= double(N);
    sd = sqrt( double(N) / double(N-1) * (sd/double(N) - mean*mean) );
}

#endif //MBL_V1_GENERIC_FUNC_H
