//
// Created by Liangsheng Zhang on 6/11/15.
//

/*
 * The function density_evol computes the evolution of a system given by density matrix.
 */

#include <iostream>
#include <iomanip>
#include "evol_op.h"
#include "evol_data.h"
#include "init_obj.h"
#include "screen_output.h"

using namespace std;

void density_evol(EvolOP* floquet, const InitObj& init_obj, EvolData& evol_data, int model = 0){

    TransitionMatrix transition;
    MatrixXcd init_density; // Initial density matrix
    vector<VectorXcd> eval; // Eigenvalues

    cout << "Diagonalize Evolution Operator." << endl;
    floquet -> Evol_Para_Init();
    floquet -> Evol_Construct();
    floquet -> Evol_Diag();
    floquet -> OP_Erase();
    floquet -> Eval(eval);

    cout << "Construct Transition Matrix." << endl;
    floquet -> Transition_Compute(transition, "Basic_Full");

    floquet -> Eigen_Erase();

    if (evol_data.evol_info.debug){
        cout << "Eigenvectors and eigenvalues:" << endl;

        for (int i=0; i < eval.size(); i++){
            cout << "Sector " << i <<" :" << endl;
            cout << "Eigenvalues:" << endl;
            complex_matrix_write(eval[i]);
            cout << endl;
        }
    }

    if (evol_data.evol_info.debug){
        cout << "Full to Basic transtion matrix:" << endl;
        complex_matrix_write(transition.Matrix("Basic_Full"));
        cout << endl;
    }

    for (int n=0; n<evol_data.evol_info.num_realization; n++){
        cout << endl;
        cout << n << "th realization:" << endl;

        InitEvolInfo init_evol_info;
        InitEvolData init_evol_data;

        cout << "Construct Initial State." << endl;
        init_density = MatrixXcd::Zero(floquet -> Get_Dim(), floquet -> Get_Dim());
        init_obj.Init_Func_C(transition, init_density, init_evol_data);

        init_evol_info.model = model;
        init_evol_info.realization = n;
        init_evol_info.debug = evol_data.evol_info.debug;

        evol_data.Init_Evol_Data(init_evol_data, init_evol_info);

        cout << "Time evolution starts." << endl;

        #pragma omp parallel num_threads(evol_data.evol_info.threads_N)
        {
            #pragma omp for
            for (int t=0; t < evol_data.evol_info.time_step; t++){
                double power = t*evol_data.evol_info.jump;

                // If time changes logarithmically
                if (evol_data.evol_info.log_time){
                    power = pow(evol_data.evol_info.log_time_jump,t);
                }

                if (power < 0){
                    cout << "Overflow happens in time evolution." << endl;
                    abort();
                }

                // Current density matrix in evec basis
                MatrixXcd density_evec(init_density.rows(), init_density.cols());
                // Current density matrix in binary basis
                VectorXcd density_basic(init_density.rows(), init_density.cols());

                vector< complex<double> > eval_power(floquet -> Get_Dim()); // power of eigenvalues

                // Here assumes the eigenvalues are eigenvalues from the unitary time evolution matrix
                int index = 0;
                for (int j=0; j < eval.size(); j++){
                    for (int k=0; k < eval[j].rows(); k++){
                        eval_power[index] = eval[j](k);
                        eval_power[index] /= abs(eval_power[index]);
                        eval_power[index] = pow(eval_power[index], power);
                        index ++;
                    }
                }


                for (int j=0; j < init_density.rows(); j++){
                    for (int k=0; k < init_density.cols(); k++){
                        density_evec(j,k) = eval_power[j] * conj(eval_power[k]) * init_density(j,k);
                    }
                }

                density_basic = transition.Matrix("Basic_Full") * density_evec
                                * transition.Matrix("Basic_Full").adjoint();

                if (density_basic.rows() != init_density.rows()){
                    cout << "Rows of current density matrix in binary basis is wrong." << endl;
                    cout << "Expected size: " << init_density.rows() << endl;
                    cout << "Obtained size: " << density_basic.rows() << endl;
                    abort();
                }

                if (density_basic.cols() != init_density.cols()){
                    cout << "Cols of current density matrix in binary basis is wrong." << endl;
                    cout << "Expected size: " << init_density.cols() << endl;
                    cout << "Obtained size: " << density_basic.cols() << endl;
                    abort();
                }

                if (evol_data.evol_info.debug){
                    cout << "Time step:" << t << endl;
                    cout << "Density Matrix in binary basis:"<<endl;
                    complex_matrix_write(density_basic);
                    cout << endl;

                    cout << "Density Matrix in evec basis:" << endl;
                    complex_matrix_write(density_evec);
                    cout << endl;
                }

                StepInfo info;
                info.realization = n;
                info.time = t;
                info.debug = evol_data.evol_info.debug;
                info.left_size = evol_data.evol_info.left_size;
                info.model = model;

                info.basis_type = "Binary";
                evol_data.Data_Compute(density_basic, info);

                info.basis_type = "Evec";
                evol_data.Data_Compute(density_evec, info);
            }
        }

        cout << "Time evolution ends." << endl;

    }

    transition.Erase_All();
}

