//
// Created by Liangsheng Zhang on 4/15/15.
//

/*
 * Implementation of functions in generic_func.h without templates to avoid duplicates when
 * included in different files.
 */

#include "generic_func.h"

using namespace std;

/*
 * Generic conjugate function
 */
double generic_conj(double val){return val;}
int generic_conj(int val){return val;}

/*
 * Generic norm function
 */
double generic_norm(double val){return abs(val);}
int generic_norm(int val){
    if (val<0) return -val;
    else return val;
}
