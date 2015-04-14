//
// Created by Liangsheng Zhang on 4/14/15.
//

#include <iostream>
#include <Eigen/Core>
#include <complex>
#include "matrix_algebra.h"

using namespace std;
using namespace Eigen;

int main(){
    MatrixXd A(2,2);
    MatrixXcd B(2,2);

    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            A(j,i) = i + j;
            B(j,i) = complex<double>(i,j);
        }
    }

    cout << "A:" << endl;
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            cout << A(i,j) << "  ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "B:" << endl;
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            cout << B(i,j) << "  ";
        }
        cout << endl;
    }
    cout << endl;

    MatrixXcd C;
    Matrix_Add(A, B, C);

    cout << "Add C:" << endl;
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            cout << C(i,j) << "  ";
        }
        cout << endl;
    }
    cout << endl;

    Matrix_Mul(A, B, C);

    cout << "Mul C:" << endl;
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            cout << C(i,j) << "  ";
        }
        cout << endl;
    }
    cout << endl;
}

