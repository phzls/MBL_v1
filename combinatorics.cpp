//
// Created by Liangsheng Zhang on 9/17/15.
//

#include <iostream>
#include <cstdlib>
#include "combinatorics.h"

using namespace std;

/*
 * Find the representation of combinotoric numbers which chooses m numbers
 * between min and max. Every number must be non-negative.
 *
 * It is done by recursively choosing the largest number.
 */
void comb_num(int min, int max, int m, vector<int>& comb, int offset ) {

    if(offset < 0){
        cout << "Offset must be non-negative." << endl;
        cout << "Offset: " << offset << endl;
        abort();
    }

    if(m==0){
        // In the last call of recursion, max can be -1
        comb.push_back(offset);
        return;
    }

    if(m<0) {
        cout << "Cannot pick negative number of numbers." << endl;
        cout << "NUmber of numbers to be picked: " << m << endl;
        abort();
    }

    if(min<0 || max<0) {
        cout << "Only non-negative numbers can be picked." << endl;
        cout << "Min: " << min << " Max: " << max << endl;
        abort();
    }

    // The lower bound which gives room for choosing smaller numbers
    int true_min = min + m-1;

    for(int i=true_min; i<=max; i++){ // The maximum number chosen
        int add = 1 << i; // The contribution from this number in the representation
        int new_offset = offset + add;
        comb_num(min, i-1, m-1, comb, new_offset);
    }
}

/*
 * Compute nCr using Pascal triangle
 */
int binomial_coef(int n, int r){
    if(n<0 || r<0){
        cout << "n and r must be positive for nCr." << endl;
        cout << "n: " << n << " r: " << r << endl;
        abort();
    }

    if(n==r) return 1;
    if(n<r) return 0; // After these checks, n>=1

    vector<int> pascal(n+1,0); // Pascal triangle
    pascal[0] = 1; // When n = 0

    // Using nCr = (n-1)C(r-1) + (n-1)Cr
    for(int i=1; i<=n; i++){
        int prev = 1; // The number in the previous position: (n-1)C(r-1)
        for(int j=1; j<i; j++){
            int temp = pascal[j];
            pascal[j] = prev + temp;
            prev = temp;
        }
        pascal[i] = 1; // The last number in the row
    }

    return pascal[r];
}
