//
// Created by Liangsheng Zhang on 4/14/15.
//

#include <iostream>
#include <Eigen/Core>
#include <complex>
#include <ctime>
#include "matrix_algebra.h"
#include "randomc.h"
#include "screen_output.h"
#include "evol_op.h"
#include "flo_evol_model.h"
#include "tasks_models.h"
#include "ham_evol_model.h"

using namespace std;
using namespace Eigen;

CRandomMersenne RanGen_mersenne(time(NULL));

TasksModels tasks_models; // Record all the tasks and methods

int EvolOP::model_num = 0;

int main(){

    clock_t time_begin = clock();
    EvolOP* floquet;

    floquet = new FloEvolHeisenRandomCosSzSectorShiftRealTau(4, 0.6, 0.8, 0, true);

    floquet -> Evol_Para_Init();
    floquet -> Evol_Construct();
    clock_t time_end = clock();
    cout << "Construction time: " << double(time_end - time_begin) / CLOCKS_PER_SEC << "s" << endl;

    time_begin = clock();
    floquet -> Evol_Diag();

    vector<MatrixXd> evec;
    vector<VectorXcd> eval;

    floquet -> Evec(evec);
    floquet -> Eval(eval);

    floquet -> Eigen_Erase();

    time_end = clock();

    cout << "Floquet time: " << double(time_end - time_begin) / CLOCKS_PER_SEC << "s" << endl;

    cout << "Eigenvalues:" << endl;
    complex_matrix_write(eval[0]);
    cout << endl;

    cout << "Eigenvectors:" << endl;
    real_matrix_write(evec[0]);
    cout << endl;

    floquet -> Eigen_Erase();

    cout << "After erasing:" << endl;
    cout << "Eigenvalues:" << endl;
    complex_matrix_write(eval[0]);
    cout << endl;

    cout << "Eigenvectors:" << endl;
    real_matrix_write(evec[0]);
    cout << endl;

    /*cout << "Eigenvalues:" << endl;
    for (int i=0; i<eval.size();i++){
        cout << "Sec " << i << endl;
        complex_matrix_write(eval[i]);
        cout << endl;
    }

    cout << "Eigenvectors:" << endl;
    for (int i=0; i<evec.size();i++){
        cout << "Sec " << i << endl;
        complex_matrix_write(evec[i]);
        cout << endl;
    }*/

    delete floquet;
    floquet = NULL;

    return 0;
}

