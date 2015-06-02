//
// Created by Liangsheng Zhang on 6/1/15.
//

#include <iostream>
#include <string>
#include "evol_op.h"
#include "parameters.h"
#include "init_obj.h"
#include "evol_data.h"
#include "tasks_models.h"

using namespace std;

/**
 ** This file implements the time evolution of multiple realizations of floquet system. There will
 ** be only one initial state in any case. If the operator has multi-sector, then the states of different
 ** sectors stored in the initial state in an order consistent with the order of sectors in eigen vector
 ** of evol_class.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.
void state_evol(EvolOP*, const InitObj&, EvolData&, int n = 0); // Evolve using state vectors. The last integer is model
// number, whichby default is 0

void multi_model_time_evolution(const AllPara& parameters){

    const string model = parameters.generic.model;
    const bool debug = parameters.generic.debug; // Whether print debug information
    const string init_func_name = parameters.evolution.init_func_name;
    const int model_num = parameters.evolution.model_num;

    EvolData evol_data(parameters);

    EvolOP* floquet;
    InitObj init_obj;

    for(int n=0; n<model_num; n++){
        cout << endl;
        cout << n << "th model" << endl;
        init_obj.init_info.size = parameters.generic.size;
        init_obj.init_info.norm_delta = 1.0e-15;
        init_obj.init_info.debug = debug;
        init_obj.init_info.init_func_name = init_func_name;

        cout << "Initialize Model." << endl;
        tasks_models.Model(model, parameters, floquet);
        init_obj.init_info.dim = floquet -> Get_Dim();

        state_evol(floquet, init_obj, evol_data, n);
    }

    cout << "Output data." << endl;

    string init_string = init_func_name;
    replace(init_string.begin(), init_string.end(),' ','_');

    string task_string = parameters.generic.task;
    replace(task_string.begin(), task_string.end(),' ','_');

    stringstream file_name;
    file_name << floquet -> Repr() << ",Task_" << task_string << ",model_num_" << model_num << ",Init_" << init_string;

    evol_data.Data_Output(parameters, file_name.str() );

    cout << endl;
    cout << endl;

    delete floquet;
    floquet = NULL;

}

