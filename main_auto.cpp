//
// Created by Liangsheng Zhang on 4/18/15.
//

#include <iostream>
#include <ctime>
#include <string>
#include "constants.h"
#include "parameters.h"
#include "tasks_models.h"
#include "randomc.h"
#include "para_model_task.h"

using namespace std;

CRandomMersenne RanGen_mersenne(time(NULL));

TasksModels tasks_models; // Record all the tasks and methods

int EvolOP::model_num = 0;

int main(int argc, char* argv[]){
    AllPara parameters;


    if (argc != 2){
        cout << "The program takes one argument: count." << endl;
        cout << "Number of argument: " << argc << endl;
        abort();
    }

    string count = argv[1];


    cout << "Count is " << count << endl;

    typedef void (*para)(AllPara&, string); // All parameters-reading functions

    map<string, para> model_para_read_map;
    map<string, para> task_para_read_map;

    model_para_read_map["Ising_Random_Simp_Flo"] = Ising_Random_Simp_Flo_Para;
    model_para_read_map["Ising_Random_Simp_Shift_Real_Flo"] = Ising_Random_Simp_Shift_Real_Flo_Para;
    model_para_read_map["Ising_Quasi_Simp_Shift_Real_Flo"] = Ising_Quasi_Simp_Shift_Real_Flo_Para;
    model_para_read_map["Ising_All_Random_Simp_Shift_Real_Flo"] = Ising_All_Random_Simp_Shift_Real_Flo_Para;

    task_para_read_map["Disorder_Transition"] = Disorder_Transition_Para;

    cout << "Generic parameters:"<< endl;
    // Obtain generic parameters
    Generic_Para(parameters, count);

    cout << "Output parameters:" << endl;
    // Obtain output parameters
    Output_Para(parameters, count);

    cout << "Model parameters:" << endl;
    // Obtain model parameters
    model_para_read_map[parameters.generic.model](parameters,count);

    cout << "Task parameters:" << endl;
    // Obtain task parameters
    task_para_read_map[parameters.generic.task](parameters, count);

    Eigen::initParallel();

    cout << "Start running." << endl;

    tasks_models.Task(parameters.generic.task)(parameters);

    cout << "Calculation finished." << endl;

    return 0;
}


