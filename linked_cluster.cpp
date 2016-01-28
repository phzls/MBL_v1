//
// Created by Liangsheng Zhang on 1/27/16.
//

#include <iostream>
#include "parameters.h"
#include "tasks_models.h"
#include "linked_cluster_class.h"
#include "methods.h"

using namespace std;

/*
 * This file implements functions related to studying coefficients of linked clusters.
 */

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

void linked_cluster(const AllPara& parameters){
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const double J_max = parameters.floquet.J_max; // Maximum J
    const double J_min = parameters.floquet.J_min; // Minimum J
    const int order_N = parameters.linked_cluster_para.order; // Number of orders
    // Total Spin Z; for now only used for Hamiltonian
    const int total_spin_z = parameters.floquet.total_spin_z;
    const int num_realization = parameters.generic.num_realizations;
    const int threads_N = parameters.generic.threads_N;
    const string model_name = parameters.generic.model;
    const bool debug = parameters.generic.debug;
    const bool time = parameters.generic.time;
    string name;

    LinkedCluster linked_cluster(parameters);

    Eigen::initParallel();

    cout << "Linked Cluster:" << endl;
    for (int i=0; i<J_N; i++){
        double J;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        cout << "J: " << J << endl;
        cout << "Start computation." << endl;

        clock_t time_begin = clock();

        ClusterData cluster_data;

        // Generating parameters for all realizations
        // For now assume all parameters are real
        cluster_data.model_para.resize(num_realization);
        for(int j=0; j<num_realization;j++){
            model_para_generation(parameters, cluster_data.model_para[j], 2*(order_N-1));
        }

        #pragma omp parallel num_threads(threads_N)
        {
            #pragma omp for
            for (int k=0; k < num_realization; k++){

                ClusterLocalInfo local_info;
                local_info.J_index = i;
                local_info.run_index = k;
                local_info.J = J;

                linked_cluster.Compute(parameters, cluster_data, local_info);
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

    linked_cluster.Output(parameters, name);
    clock_t out_end = clock();

    if (time){
        cout << "Output time: " << double(out_end - out_begin) / CLOCKS_PER_SEC << "s" << endl;
        cout << endl;
    }
}





