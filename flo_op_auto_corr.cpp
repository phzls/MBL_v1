//
// Created by Liangsheng Zhang on 9/24/15.
//

/**
 ** This program computes autocorrelation of some time independent operator for
 ** Floquet system. Particularly, let the operator be O, and U be the time evolution
 ** operator of the system, then it computes Tr{ O*(U^{dagger}^n)*O*(U^n) } for certain
 ** n. It will also do average over random ensembles.
 **
 ** The O is constructed beforehand, and converts in the basis of U. The trace is
 ** computed by a double sum without invoking any matrix multiplication.
 **/

#include <iostream>
#include <ctime>
#include "parameters.h"
#include "tasks_models.h"
#include "screen_output.h"
#include "evol_op.h"
#include "flo_op_auto_corr.h"

using namespace std;

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void flo_op_auto_corr(const AllPara& parameters){
    const int para_pts = parameters.flo_op_auto_corr.para_pts; // Number of points for parameter
    const int num_realizations = parameters.generic.num_realizations;
    const int time_pts = parameters.flo_op_auto_corr.time_pts;
    const int threads_N = parameters.generic.threads_N;
    const string model_name = parameters.generic.model;
    const bool debug = parameters.generic.debug;
    const bool time = parameters.generic.time;
    string name;

    OpAutoCorr op_auto_corr(parameters);
    OpCorrLocalPara op_corr_local_para(parameters);

    AllPara local_parameters(parameters);

    string para_name = op_corr_local_para.Para_Name();

    for(int i=0; i<para_pts;i++){
        double para = op_corr_local_para.Para_Update(i, local_parameters);

        cout << para_name << ": " << para << endl;
        cout << "Start computation." << endl;

        clock_t time_begin = clock();

        for(int j=0; j<num_realizations; j++){
            EvolOP* model;
            tasks_models.Model(model_name, local_parameters, model);
            model -> Evol_Para_Init();
            model -> Evol_Construct();
            model -> Evol_Diag(true); // Eigenvectors are kept
            model -> Eval(op_auto_corr.eval);

            if (i == 0 && j == 0) {
                name = model -> Type();
            }

            op_auto_corr.SetUp(local_parameters, model); // Set up the operators

            model -> OP_Erase();
            model -> Eigen_Erase();

            delete model;
            model = NULL;

            if (debug){
                cout << "Realization " << j << endl;

                cout << "Eigenvalues:" << endl;

                for (int l=0; l<op_auto_corr.eval.size(); l++){
                    cout << "Sector " << l << endl;
                    complex_matrix_write(op_auto_corr.eval[l]);
                    cout << endl;
                }
            }

            #pragma omp parallel num_threads(threads_N)
            {
                #pragma omp for
                for(int k=0; k<time_pts; k++){
                    OpCorrLocalInfo op_corr_local_info;
                    op_corr_local_info.para_num = i;
                    op_corr_local_info.realization_num = j;
                    op_corr_local_info.time_step = k;

                    op_auto_corr.Compute(parameters, op_corr_local_info);
                }
            }
        }

        clock_t time_end = clock();

        if (time){
            cout << "Computation time: " << double(time_end - time_begin) / CLOCKS_PER_SEC << "s" << endl;
            cout << endl;
        }
    }

    cout << "Output data." << endl;
    clock_t out_begin = clock();

    // Add J to floquet name
    if( name.find("Flo") != string::npos ){
        stringstream new_name;
        new_name << name << ",W=" << parameters.floquet.J;
        name = new_name.str();
    }

    op_auto_corr.Output(parameters, op_corr_local_para, name);
    clock_t out_end = clock();

    if (time){
        cout << "Output time: " << double(out_end - out_begin) / CLOCKS_PER_SEC << "s" << endl;
        cout << endl;
    }

}




