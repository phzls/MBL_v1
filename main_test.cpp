//
// Created by Liangsheng Zhang on 4/16/15.
//

#include <iostream>
#include <ctime>
#include <string>
#include "constants.h"
#include "parameters.h"
#include "tasks_models.h"
#include "randomc.h"

using namespace std;

CRandomMersenne RanGen_mersenne(time(NULL));

TasksModels tasks_models; // Record all the tasks and methods

int EvolOP::model_num = 0;

int main() {
    AllPara parameters;




//===================================================================================





    parameters.generic.task = "Multi_Model_Time_Evolution";
    parameters.generic.model = "Ising_Random_Simp_Shift_Real_Flo";

    parameters.generic.size = 4; // System size
    parameters.generic.num_realizations = 1; // Number of realizations
    parameters.generic.threads_N = 1; // Number of threads in openmp
    parameters.generic.debug = true; // Whether output debug information
    parameters.generic.version = 1; // Version of the output
    parameters.generic.time = true; // Whether the program is timed

    parameters.output.width = 30; // Width for spacing in output files
    parameters.output.filename_output = true; // Whether print out file names

    parameters.floquet.J_N = 30; // Number of points of coupling strength
    parameters.floquet.J_min = 0.1; // Minimum J
    parameters.floquet.J_max = 0.9; // Maximum J
    parameters.floquet.tau = 0.8; // Time step size
    parameters.floquet.J = 0.6;

    parameters.markov.K = 0.8; // Coupling strength to the bath in markov models




//===============================================================================




    parameters.evolution.time_step = 10; // Number of time steps
    parameters.evolution.step_size = parameters.floquet.tau; // Time step size
    parameters.evolution.init_func_name = "Left_Spin_Random"; // Initial state name

    parameters.evolution.evol_compute["Entropy_Per_Model"] = false;
    parameters.evolution.evol_compute["Leftmost_Spin_Per_Model"] = true;

    parameters.evolution.sample_detail = true; // Determine whether print out all sample values

    parameters.evolution.model_num = 1; // Number of models for evolution
    // If partition the chain to two halves, the size of left part
    parameters.evolution.left_size = parameters.generic.size / 2;
    parameters.evolution.jump = 1; // jump of time points in evolution

    parameters.evolution.markov_time_jump = 10; // Time jump in the markov process for flipping
    // True time: time_step * markov_time_jump
    parameters.evolution.markov_jump = true; // Determine whether there will be markov_time_jump
    if (parameters.generic.task.find("Markov") == string::npos &&
        parameters.generic.task.find("Single Model") == string::npos) {
        // The task has nothing to do with markov process
        parameters.evolution.markov_jump = false;
    }

    parameters.evolution.log_time = true; // whehter time changes logarithmically
    parameters.evolution.log_time_jump = 2; // The base for time change logarithmically

    // The number gives the index of leftmost spin z value
    parameters.evolution.leftmost_spin_z_index = 3;

    // Multiple sets of initial conditions
    int leftmost_spin_z_index_set[] = {1, 2, 3, 4, 5, 63, 64};
    parameters.multi_ini_para.leftmost_spin_z_index_set.assign(leftmost_spin_z_index_set,
                                                               leftmost_spin_z_index_set +
                                                               sizeof(leftmost_spin_z_index_set) / sizeof(int));

    // Whether output the evolution of leftmot spin z for all eigenstates in full_leftmost_spin_z
    parameters.multi_ini_para.easy_full_leftmost_spin_z = false;

    // Threshold in time evolution of leftmost spin z for non-zero values
    parameters.multi_ini_para.non_zero_threshold = 0.001;

    parameters.evolution.evol_way = "matrix";


//=================================================================================





    // Methods to be called under single model
    parameters.single_model.single_model_compute["Flo Chain End Sigma Z"] = false;
    parameters.single_model.single_model_compute["Flo Evolution Simple Markov"] = true;




//===================================================================================




    // Methods to be called for studying transition of floquet systems from thermal to localization
    parameters.transition.flo_transition_compute["ZZ_Correlation_Square"] = true; // End-to-end sigma_z X sigma_z
    // correlation square

    parameters.transition.flo_transition_compute["Entropy_Variance"] = true; // Entropy variance for all eigenstates

    parameters.transition.flo_transition_compute["ZZ_Time_Correlation"] = true; // End-to-end sigma_z X sigma_z
    // time correlation

    parameters.transition.flo_transition_compute["ZZ_Time_Correlation_Components"] = false; // End-to-end
    // sigma_z X sigma_z time correlation components. The first row in output is text header

    parameters.transition.flo_transition_compute["ZZ_All_Correlation-Square"] = false; // zz correlation square
    // at all distances

    parameters.transition.flo_transition_compute["Log_ZZ_Correlation_Square"] = false; // zz correlation suqare with
    // logarithm taken with each eigenstate first

    parameters.transition.flo_transition_compute["Log_ZZ_All_Correlation_Square"] = false; // log zz correlation square
    // at all distances

    parameters.transition.flo_transition_compute["ZZ_Correlation_Square_All_Sample"] = false; // zz correlation square
    // for all samples

    parameters.transition.flo_transition_compute["ZZ_Time_Correlation_All_Sample"] = false; // zz time correlation for
    // all samples

    parameters.transition.flo_transition_compute["ZZ_All_Time_Correlation"] = false; // zz time correlation at all
    // distances




//============================================================================================



    // Methods to be called for studying eigenstate properties of floquet systems
    parameters.eigenvec.flo_eigen_compute["ZZ Correlation Square"] = false; // End-to-end sigma_z sigma_z correlation
    // square
    parameters.eigenvec.flo_eigen_compute["Eigenvectors"] = true; // End-to-end sigma_z sigma_z correlation
    // square



//============================================================================================




    Eigen::initParallel();

    tasks_models.Task(parameters.generic.task)(parameters);

    cout << "Calculation finished." << endl;

    return 0;
}


