//
// Created by Liangsheng Zhang on 9/17/15.
//

#ifndef MBL_V1_COMBINATORICS_H
#define MBL_V1_COMBINATORICS_H

#include <vector>

using namespace std;

/**
 ** This file implements functions related to combinatorics
 **/


/*
 * This function generates combinatoric number in the following sense:
 * it chooses m numbers between min and max, both ends inclusive.
 * For each combination, we can think of picking up 1 in corresponding positions,
 * which corresponds to a binary representation of a number. An overal offest can
 * also be added to this representation, which is given by offset.
 * These (possibly with offset) representation numbers are stored in vector comb.
 * All numbers must be non-negative.
 */
void comb_num(int min, int max, int m, vector<int>& comb, int offset = 0);

/*
 * Compute combinatoric number nCr, where n and r must be non-negative integers
 */
int binomial_coef(int n, int r);


#endif //MBL_V1_COMBINATORICS_H
