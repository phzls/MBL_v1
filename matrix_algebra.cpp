//
// Created by Liangsheng Zhang on 4/14/15.
//

/*
 * Implement functions defined in matrix_algebra.h
 */

#include <complex>
#include "matrix_algebra.h"

using namespace std;
using namespace Eigen;

void Matrix_Add(const MatrixXd& A, const MatrixXd& B, MatrixXd& C){
    Matrix_Add_Sub_Check(A,B);
    C = A + B;
}

void Matrix_Add(const MatrixXd& A, const MatrixXd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);
    MatrixXd temp;
    temp = A + B;

    C = MatrixXcd::Zero(temp.rows(), temp.cols());

    for (int i=0; i<C.cols(); i++){
        for (int j=0; j<C.rows(); j++){
            C(j,i) = complex<double>(temp(j,i),0);
        }
    }
}

void Matrix_Add(const MatrixXd& A, const MatrixXcd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);
    MatrixXcd temp(A.cols(),A.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(A(j,i),0);
        }
    }

    C = temp + B;
}

void Matrix_Add(const MatrixXcd& A, const MatrixXd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);
    MatrixXcd temp(B.cols(),B.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(B(j,i),0);
        }
    }

    C = A + temp;
}

void Matrix_Add(const MatrixXcd& A, const MatrixXcd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);

    C = A + B;
}

void Matrix_Sub(const MatrixXd& A, const MatrixXd& B, MatrixXd& C){
    Matrix_Add_Sub_Check(A,B);
    C = A - B;
}

void Matrix_Sub(const MatrixXd& A, const MatrixXd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);
    MatrixXd temp;
    temp = A - B;

    C = MatrixXcd::Zero(temp.rows(), temp.cols());

    for (int i=0; i<C.cols(); i++){
        for (int j=0; j<C.rows(); j++){
            C(j,i) = complex<double>(temp(j,i),0);
        }
    }
}

void Matrix_Sub(const MatrixXd& A, const MatrixXcd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);
    MatrixXcd temp(A.cols(),A.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(A(j,i),0);
        }
    }

    C = temp - B;
}

void Matrix_Sub(const MatrixXcd& A, const MatrixXd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);
    MatrixXcd temp(B.cols(),B.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(B(j,i),0);
        }
    }

    C = A - temp;
}

void Matrix_Sub(const MatrixXcd& A, const MatrixXcd& B, MatrixXcd& C){
    Matrix_Add_Sub_Check(A,B);

    C = A - B;
}

void Matrix_Mul(const MatrixXd& A, const MatrixXd& B, MatrixXd& C){
    Matrix_Mul_Check(A,B);

    C = A * B;
}

void Matrix_Mul(const MatrixXd& A, const MatrixXd& B, MatrixXcd& C){
    Matrix_Mul_Check(A,B);
    MatrixXd temp;
    temp = A * B;

    C = MatrixXcd::Zero(temp.rows(), temp.cols());

    for (int i=0; i<C.cols(); i++){
        for (int j=0; j<C.rows(); j++){
            C(j,i) = complex<double>(temp(j,i),0);
        }
    }
}

void Matrix_Mul(const MatrixXd& A, const MatrixXcd& B, MatrixXcd& C){
    Matrix_Mul_Check(A,B);
    MatrixXcd temp(A.cols(),A.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(A(j,i),0);
        }
    }

    C = temp * B;
}

void Matrix_Mul(const MatrixXcd& A, const MatrixXd& B, MatrixXcd& C){
    Matrix_Mul_Check(A,B);
    MatrixXcd temp(B.cols(),B.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(B(j,i),0);
        }
    }

    C = A * temp;
}

void Matrix_Mul(const MatrixXcd& A, const MatrixXcd& B, MatrixXcd& C){
    Matrix_Mul_Check(A,B);

    C = A * B;
}

void Matrix_Add(const MatrixXd& A, const MatrixXd& B, VectorXd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    C = A + B;
}

void Matrix_Add(const MatrixXd& A, const MatrixXd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    MatrixXd temp;
    temp = A + B;

    C = MatrixXcd::Zero(C.rows(), C.cols());

    for (int i=0; i<C.cols(); i++){
        for (int j=0; j<C.rows(); j++){
            C(j,i) = complex<double>(temp(j,i),0);
        }
    }
}

void Matrix_Add(const MatrixXd& A, const MatrixXcd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    MatrixXcd temp(A.cols(),A.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(A(j,i),0);
        }
    }

    C = temp + B;
}

void Matrix_Add(const MatrixXcd& A, const MatrixXd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    MatrixXcd temp(B.cols(),B.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(B(j,i),0);
        }
    }

    C = A + temp;
}

void Matrix_Add(const MatrixXcd& A, const MatrixXcd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    C = A + B;
}

void Matrix_Sub(const MatrixXd& A, const MatrixXd& B, VectorXd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    C = A - B;
}

void Matrix_Sub(const MatrixXd& A, const MatrixXd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    MatrixXd temp;
    temp = A - B;

    C = MatrixXcd::Zero(C.rows(), C.cols());

    for (int i=0; i<C.cols(); i++){
        for (int j=0; j<C.rows(); j++){
            C(j,i) = complex<double>(temp(j,i),0);
        }
    }
}

void Matrix_Sub(const MatrixXd& A, const MatrixXcd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    MatrixXcd temp(A.cols(),A.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(A(j,i),0);
        }
    }

    C = temp - B;
}

void Matrix_Sub(const MatrixXcd& A, const MatrixXd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    MatrixXcd temp(B.cols(),B.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(B(j,i),0);
        }
    }

    C = A - temp;
}

void Matrix_Sub(const MatrixXcd& A, const MatrixXcd& B, VectorXcd& C){
    Matrix_Add_Sub_Check(A,B);
    Vector_Check(A);
    Vector_Check(B);

    C = A - B;
}

void Matrix_Mul(const MatrixXd& A, const MatrixXd& B, VectorXd& C){
    Matrix_Mul_Check(A,B);
    Mul_Vector_Check(A,B);

    C = A * B;
}

void Matrix_Mul(const MatrixXd& A, const MatrixXd& B, VectorXcd& C){
    Matrix_Mul_Check(A,B);
    Mul_Vector_Check(A,B);

    MatrixXd temp;
    temp = A * B;

    C = MatrixXcd::Zero(C.rows(), C.cols());

    for (int i=0; i<C.cols(); i++){
        for (int j=0; j<C.rows(); j++){
            C(j,i) = complex<double>(temp(j,i),0);
        }
    }
}

void Matrix_Mul(const MatrixXd& A, const MatrixXcd& B, VectorXcd& C){
    Matrix_Mul_Check(A,B);
    Mul_Vector_Check(A,B);

    MatrixXcd temp(A.cols(),A.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(A(j,i),0);
        }
    }

    C = temp * B;
}

void Matrix_Mul(const MatrixXcd& A, const MatrixXd& B, VectorXcd& C){
    Matrix_Mul_Check(A,B);
    Mul_Vector_Check(A,B);

    MatrixXcd temp(B.cols(),B.rows());

    for (int i=0; i<temp.cols(); i++){
        for (int j=0; j<temp.rows(); j++){
            temp(j,i) = complex<double>(B(j,i),0);
        }
    }

    C = A * temp;
}

void Matrix_Mul(const MatrixXcd& A, const MatrixXcd& B, VectorXcd& C){
    Matrix_Mul_Check(A,B);
    Mul_Vector_Check(A,B);

    C = A * B;
}