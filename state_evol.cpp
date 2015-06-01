//
// Created by Liangsheng Zhang on 6/1/15.
//

/*
 * The function state_evol computes the evolution of a system given by state vectors.
 */

#include <iostream>
#include "evol_op.h"
#include "evol_data.h"
#include "init_obj.h"
#include "screen_output.h"

using namespace std;

void state_evol(EvolOP* floquet, const InitObj& init_obj, EvolData& evol_data){

    TransitionMatrix transition;
    VectorXcd init_state; // Initial state
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

        cout << "Construct Initial State." << endl;
        init_state = VectorXcd::Zero(floquet -> Get_Dim());
        init_obj.Init_Func(transition, init_state);

        cout << "Time evolution starts." << endl;

        #pragma omp parallel num_threads(threads_N)
        {
            #pragma omp for
            for (int t=0; t < evol_data.evol_info.time_step; t++){
                long long int power = t*evol_data.evol_info.jump;

                // If time changes logarithmically
                if (evol_data.evol_info.log_time){
                    power = pow(evol_data.evol_info.log_time_jump,t);
                }

                if (power < 0){
                    cout << "Overflow happens in time evolution." << endl;
                    abort();
                }

                VectorXcd state_evec(init_state.size()); // Current state in evec basis
                VectorXcd state_basic(init_state.size()); // Current state in binary basis

                int index = 0;
                for (int j=0; j < eval.size(); j++){
                    for (int k=0; k < eval[j].rows(); k++){
                        complex<double> eigenvalue = eval[j](k);
                        state_evec(index) = pow(eigenvalue, power) * init_state(index);
                        index ++;
                    }
                }

                state_basic = transition.Matrix("Basic_Full") * state_evec;

                if (state_basic.size() != init_state.size()){
                    cout << "Size of current state in binary basis is wrong." << endl;
                    cout << "Expected size: " << init_state.size() << endl;
                    cout << "Obtained size: " << state_basic.size() << endl;
                    abort();
                }

                if (evol_data.evol_info.debug){
                    cout << "Time step:" << t << endl;
                    cout << "State in binary basis:"<<endl;
                    complex_matrix_write(state_basic);
                    cout << endl;

                    cout << "State in evec basis:" << endl;
                    complex_matrix_write(state_evec);
                    cout << endl;
                }

                StepInfo info;
                info.realization = n;
                info.time = t;
                info.debug = evol_data.evol_info.debug;
                info.left_size = evol_data.evol_info.left_size;

                evol_data.Data_Compute(state_basic, info);
            }
        }

        cout << "Time evolution ends." << endl;

    }

    transition.Erase_All();
}
