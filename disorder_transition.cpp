//
// Created by Liangsheng Zhang on 4/15/15.
//

#include <iostream>
#include <ctime>
#include "parameters.h"
#include "tasks_models.h"
#include "disorder_model_transition.h"
#include "screen_output.h"
#include "evol_op.h"
#include "sort_comparator.h"
#include "disorder_transition_func.h"

using namespace std;

/**
 ** This file implements functions related to studying transition of floquet models from thermal to localization.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void disorder_transition(const AllPara& parameters){
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const double J_max = parameters.floquet.J_max; // Maximum J
    const double J_min = parameters.floquet.J_min; // Minimum J
    const int total_spin_z = parameters.floquet.total_spin_z; // Total Spin Z; for now only used for Hamiltonian
    const int num_realization = parameters.generic.num_realizations;
    const int threads_N = parameters.generic.threads_N;
    const string model_name = parameters.generic.model;
    const bool debug = parameters.generic.debug;
    const bool time = parameters.generic.time;
    const bool mid_half_spectrum = parameters.transition.mid_half_spectrum; // Whether only use middle half spectrum
    string name;
    bool is_ham; // Whether the model is Hamiltonian

    DisorderModelTransition disorder_model_transition(parameters);

    AllPara local_parameters(parameters); // Local parameters which can be changed

    Eigen::initParallel();

    cout << "Disorder transition:" << endl;
    for (int i=0; i<J_N; i++){
        double J;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        local_parameters.floquet.J = J;
        local_parameters.floquet.total_spin_z = total_spin_z;

        vector<EvolOP*> models(num_realization);

        for (int k=0; k<num_realization; k++){
            tasks_models.Model(model_name, local_parameters, models[k]);
            models[k] -> Evol_Para_Init();
        }

        if (i == 0) {
            name = models[0] -> Type();
            is_ham = (name.find("Hamiltonian") != string::npos);
            if(is_ham){
                if(mid_half_spectrum) name += ",mid_half_spec";
                else name += ",full_spec";
            }
        }

        cout << "J: " << local_parameters.floquet.J << endl;
        cout << "Start computation." << endl;

        clock_t time_begin = clock();

        #pragma omp parallel num_threads(threads_N)
        {
            #pragma omp for
            for (int k=0; k < num_realization; k++){

                DisorderLocalInfo local_info;
                local_info.J_index = i;
                local_info.realization_index = k;
                local_info.is_ham = is_ham;
                local_info.mid_half_spectrum = mid_half_spectrum;

                models[k] -> Evol_Construct();

                if (disorder_model_transition.Op_Diag()) {
                    // Diagonalize the matrix
                    models[k] -> Evol_Diag(disorder_model_transition.Op_Evec_Keep());

                    if (disorder_model_transition.Eval_Real() && (models[k]->Eval_Type() == "Real")) {
                        models[k] -> Eval(local_info.eval_real);
                        local_info.eval_type_real = true;
                    }
                    else{
                        models[k] -> Eval(local_info.eval_complex);
                        local_info.eval_type_real = false;
                    }

                    if (disorder_model_transition.Op_Evec_Keep()){
                        // Keep eigenvectors
                        if (disorder_model_transition.Evec_Real() && (models[k]->Evec_Type() == "Real")) {
                            models[k] -> Evec(local_info.evec_real);
                            local_info.evec_type_real = true;
                        }
                        else{
                            models[k] -> Evec(local_info.evec_complex);
                            local_info.evec_type_real = false;
                        }
                    }

                    models[k] -> Eigen_Erase();

                    // If it's Hamiltonian, then I need to only take center half in each sector
                    if( is_ham && mid_half_spectrum ){
                        if(!local_info.eval_type_real){
                            cout << "Eigenvalues of Hamiltonian must be real" << endl;
                        }

                        if(debug) cout << "Original Eigenvalues and Eigenvectors:" << endl;
                        for(int l=0; l<local_info.eval_real.size(); l++){
                            // Second entry acts as an index
                            const int sec_size = local_info.eval_real[l].size();

                            vector< pair<double, int> > pair_index(sec_size);
                            for(int m=0; m < sec_size; m++){
                                pair_index[m] = pair<double, int>(local_info.eval_real[l](m), m);
                            }

                            if(debug){
                                cout << "Sector " << l << endl;
                                cout << "Eigenvalues:" << endl;
                                real_matrix_write(local_info.eval_real[l]);
                                cout << endl;
                            }

                            sort(pair_index.begin(), pair_index.end(), pair_first_less_comparator<double, int>);

                            // Obtain the eigenvalue above and below the middle half
                            local_info.lower_eval = pair_index[sec_size/4-1].first;
                            local_info.upper_eval = pair_index[sec_size-sec_size/4].first;

                            // Take only the middle half eigenvalues
                            Middle_Half_Extract(pair_index, local_info.eval_real[l]);

                            // Take only the middle half eigenvectors if they are computed
                            if(disorder_model_transition.Op_Evec_Keep()){
                                if(debug) cout << "Eigenvectors:" << endl;
                                if(local_info.evec_type_real){
                                    if(debug){
                                        real_matrix_write(local_info.evec_real[l]);
                                        cout << endl;
                                    }

                                    Middle_Half_Extract(pair_index, local_info.evec_real[l]);
                                }
                                else{
                                    if(debug){
                                        complex_matrix_write(local_info.evec_complex[l]);
                                        cout << endl;
                                    }

                                    Middle_Half_Extract(pair_index, local_info.evec_complex[l]);
                                }
                            }
                        }
                    }
                }

                if (!disorder_model_transition.Op_Keep()) models[k] -> OP_Erase();

                if (debug){
                    cout << "Realization " << k << endl;

                    if (disorder_model_transition.Op_Diag()){
                        cout << "Eigenvalues:" << endl;

                        if (disorder_model_transition.Eval_Real() && (models[k]->Eval_Type() == "Real")){
                            for (int l=0; l<local_info.eval_real.size(); l++){
                                cout << "Sector " << l << endl;
                                real_matrix_write(local_info.eval_real[l]);
                                cout << endl;
                            }
                        }
                        else{
                            for (int l=0; l<local_info.eval_complex.size(); l++){
                                cout << "Sector " << l << endl;
                                complex_matrix_write(local_info.eval_complex[l]);
                                cout << endl;
                            }
                        }

                        if(is_ham){
                            cout << "Lower eval: " << local_info.lower_eval << endl;
                            cout << "Upper eval: " << local_info.upper_eval << endl;
                            cout << endl;
                        }

                        if (disorder_model_transition.Op_Evec_Keep()){
                            cout << "Eigenvectors:" << endl;

                            if (disorder_model_transition.Evec_Real() && (models[k]->Evec_Type() == "Real")){
                                for (int l=0; l<local_info.evec_real.size(); l++){
                                    cout << "Sector " << l << endl;
                                    real_matrix_write(local_info.evec_real[l]);
                                    cout << endl;
                                }
                            }
                            else{
                                for (int l=0; l<local_info.evec_complex.size(); l++){
                                    cout << "Sector " << l << endl;
                                    complex_matrix_write(local_info.evec_complex[l]);
                                    cout << endl;
                                }
                            }
                        }
                    }
                }

                disorder_model_transition.Compute(parameters, models[k], local_info);

                delete models[k];
                models[k] = NULL;
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

    // Add tau value if it can be chosen
    if( name.find("Tau") != string::npos ){
        stringstream new_name;
        new_name << name << ",tau=" << parameters.floquet.tau;
        name = new_name.str();
    }

    // Add alpha value for continuous modified model
    if( name.find("Con_") != string::npos ){
        stringstream new_name;
        new_name << name << ",alpha=" << parameters.floquet.alpha;
        name = new_name.str();
    }

    disorder_model_transition.Output(parameters, name);
    clock_t out_end = clock();

    if (time){
        cout << "Output time: " << double(out_end - out_begin) / CLOCKS_PER_SEC << "s" << endl;
        cout << endl;
    }
}



