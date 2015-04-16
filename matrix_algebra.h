//
// Created by Liangsheng Zhang on 4/14/15.
//

#ifndef MATRIX_ALGEBRA_H
#define MATRIX_ALGEBRA_H

/*
 * This file defines functions which handle general matrix algebra based on Eigen.
 */

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void Matrix_Add(const MatrixXd&, const MatrixXd&, MatrixXd&);
void Matrix_Add(const MatrixXd&, const MatrixXd&, MatrixXcd&);
void Matrix_Add(const MatrixXd&, const MatrixXcd&, MatrixXcd&);
void Matrix_Add(const MatrixXcd&, const MatrixXd&, MatrixXcd&);
void Matrix_Add(const MatrixXcd&, const MatrixXcd&, MatrixXcd&);

void Matrix_Add(const MatrixXd&, const MatrixXd&, VectorXd&);
void Matrix_Add(const MatrixXd&, const MatrixXd&, VectorXcd&);
void Matrix_Add(const MatrixXd&, const MatrixXcd&, VectorXcd&);
void Matrix_Add(const MatrixXcd&, const MatrixXd&, VectorXcd&);
void Matrix_Add(const MatrixXcd&, const MatrixXcd&, VectorXcd&);

void Matrix_Sub(const MatrixXd&, const MatrixXd&, MatrixXd&);
void Matrix_Sub(const MatrixXd&, const MatrixXd&, MatrixXcd&);
void Matrix_Sub(const MatrixXd&, const MatrixXcd&, MatrixXcd&);
void Matrix_Sub(const MatrixXcd&, const MatrixXd&, MatrixXcd&);
void Matrix_Sub(const MatrixXcd&, const MatrixXcd&, MatrixXcd&);

void Matrix_Sub(const MatrixXd&, const MatrixXd&, VectorXd&);
void Matrix_Sub(const MatrixXd&, const MatrixXd&, VectorXcd&);
void Matrix_Sub(const MatrixXd&, const MatrixXcd&, VectorXcd&);
void Matrix_Sub(const MatrixXcd&, const MatrixXd&, VectorXcd&);
void Matrix_Sub(const MatrixXcd&, const MatrixXcd&, VectorXcd&);

void Matrix_Mul(const MatrixXd&, const MatrixXd&, MatrixXd&);
void Matrix_Mul(const MatrixXd&, const MatrixXd&, MatrixXcd&);
void Matrix_Mul(const MatrixXd&, const MatrixXcd&, MatrixXcd&);
void Matrix_Mul(const MatrixXcd&, const MatrixXd&, MatrixXcd&);
void Matrix_Mul(const MatrixXcd&, const MatrixXcd&, MatrixXcd&);

void Matrix_Mul(const MatrixXd&, const MatrixXd&, VectorXd&);
void Matrix_Mul(const MatrixXd&, const MatrixXd&, VectorXcd&);
void Matrix_Mul(const MatrixXd&, const MatrixXcd&, VectorXcd&);
void Matrix_Mul(const MatrixXcd&, const MatrixXd&, VectorXcd&);
void Matrix_Mul(const MatrixXcd&, const MatrixXcd&, VectorXcd&);

template <class T1, class T2>
void Matrix_Mul_Check(const MatrixBase<T1>& A, const MatrixBase<T2>& B){
    // Check dimensionality for multiplication
    if (A.cols() != B.rows()){
        cout << "Dimension does not match for matrix multiplication." << endl;
        cout << "First matrix col num: " << A.cols() << endl;
        cout << "Second matrix row num: " << B.rows() << endl;
        abort();
    }
}

template <class T1, class T2>
void Matrix_Add_Sub_Check(const MatrixBase<T1>& A, const MatrixBase<T2>& B){
    // Check dimensionality for addition and subtraction
    if (A.cols() != B.cols()){
        cout << "Col dimension does not match for matrix multiplication." << endl;
        cout << "First matrix col num: " << A.cols() << endl;
        cout << "Second matrix col num: " << B.cols() << endl;
        abort();
    }

    if (A.rows() != B.rows()){
        cout << "Row dimension does not match for matrix multiplication." << endl;
        cout << "First matrix row num: " << A.rows() << endl;
        cout << "Second matrix row num: " << B.rows() << endl;
        abort();
    }
}

template <class T>
void Vector_Check(const MatrixBase<T>& A){
    // Check if it is a vector
    if (A.cols() != 1 && A.rows() != 1){
        cout << "Matrix is not of vector shape." << endl;
        cout << "Row num: " << A.rows() << endl;
        cout << "Col num: " << A.cols() << endl;
        abort();
    }
}

template <class T1, class T2>
void Mul_Vector_Check(const MatrixBase<T1>& A,  const MatrixBase<T2>& B){
    // Check if the multiplication results in a vector
    if (A.rows() != 1 && B.cols() != 1){
        cout << "Multiplication does not give a vector." << endl;
        cout << "First matrix row num: " << A.rows() << endl;
        cout << "Second matrix col num: " << B.cols() << endl;
        abort();
    }
}

#endif //MATRIX_ALGEBRA_H
