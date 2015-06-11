//
// Created by Liangsheng Zhang on 6/1/15.
//

#include <iostream>
#include "evol_op.h"
#include "parameters.h"
#include "init_obj.h"
#include "evol_data.h"
#include "tasks_models.h"

using namespace std;

/**
 ** This file implements the time evolution of a floquet system. There will be only one initial
 ** state in any case. If the operator has multi-sector, then the states of different sectors
 ** stored in the initial state in an order consistent with the order of sectors in eigen vector
 ** of evol_class.
 **/

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.
void state_evol(EvolOP*, const InitObj&, EvolData&, int n=0); // Evolve using state vectors. The last integer is model
// number, whichby default is 0
void density_evol(EvolOP*, const InitObj&, EvolData&, int n = 0); // Evolve using density matrix. The last integer is model
// number, which by default is 0

void single_model_time_evolution(const AllPara& parameters){

    const string model = parameters.generic.model;
    const bool debug = parameters.generic.debug; // Whether print debug information
    const string init_func_name = parameters.evolution.init_func_name;
    const string evol_way = parameters.evolution.evol_way;

    if (parameters.evolution.model_num != 1){
        cout << "Single model time evolution should only have one model realization." << endl;
        cout << "Obtained number of models:" << parameters.evolution.model_num << endl;
        abort();
    }

    EvolOP* floquet;
    EvolData evol_data(parameters);
    InitObj init_obj;

    init_obj.init_info.size = parameters.generic.size;
    init_obj.init_info.norm_delta = 1.0e-15;
    init_obj.init_info.debug = debug;
    init_obj.init_info.init_func_name = init_func_name;

    cout << "Initialize Model." << endl;
    tasks_models.Model(model, parameters, floquet);
    init_obj.init_info.dim = floquet -> Get_Dim();

    if (evol_way == "vector") state_evol(floquet, init_obj, evol_data);
    else if (evol_way == "matrix") density_evol(floquet, init_obj, evol_data);
    else{
        cout << "Evol way " << evol_way << " is not recognized." << endl;
        abort();
    }

    cout << "Output data." << endl;

    string init_string = init_func_name;
    replace(init_string.begin(), init_string.end(),' ','_');

    string task_string = parameters.generic.task;
    replace(task_string.begin(), task_string.end(),' ','_');

    evol_data.Data_Output(parameters, floquet -> Repr() + ",Task_" + task_string + ",Init_" + init_string );

    cout << endl;
    cout << endl;

    delete floquet;
    floquet = NULL;

}

