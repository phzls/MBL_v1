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
    model_para_read_map["Ising_All_Quasi_Simp_Shift_Real_Flo"] = Ising_All_Quasi_Simp_Shift_Real_Flo_Para;
    model_para_read_map["Ising_Random_Simp_Shift_Cos_Real_Flo"] = Ising_Random_Simp_Shift_Cos_Real_Flo_Para;
    model_para_read_map["Ising_Random_Simp_Shift_Cos_Real_Tau_Flo"] = Ising_Random_Simp_Shift_Cos_Real_Tau_Flo_Para;
    model_para_read_map["Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Flo"] = Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Flo_Para;
    model_para_read_map["Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Flo"] = Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Flo_Para;
    model_para_read_map["Heisen_Random_Cos_Sz_Sector_Modified_Tau_Flo"] = Heisen_Random_Cos_Sz_Sector_Modified_Tau_Flo_Para;
    model_para_read_map["Heisen_Quasi_Sz_Sector_Modified_Tau_Flo"] = Heisen_Quasi_Sz_Sector_Modified_Tau_Flo_Para;
    model_para_read_map["XXZ_Gaussian_Random_Shift_Real_Flo"] = XXZ_Gaussian_Random_Shift_Real_Flo_Para;
    model_para_read_map["XXZ_Uniform_Random_Shift_Real_Flo"] = XXZ_Uniform_Random_Shift_Real_Flo_Para;
    model_para_read_map["XXZ_Gaussian_Z_Random_Shift_Real_Flo"] = XXZ_General_Z_Random_Shift_Real_Flo_Para;
    model_para_read_map["XXZ_Uniform_Z_Random_Shift_Real_Flo"] = XXZ_General_Z_Random_Shift_Real_Flo_Para;
    model_para_read_map["XXZ_Gaussian_Random_Field_Shift_Real_Flo"] = XXZ_General_Random_Field_Shift_Real_Flo_Para;
    model_para_read_map["XXZ_Uniform_Random_Field_Shift_Real_Flo"] = XXZ_General_Random_Field_Shift_Real_Flo_Para;

    model_para_read_map["Heisen_Random_Cos_Sz_Sector_Ham"] = Heisen_Random_Cos_Sz_Sector_Ham_Para;
    model_para_read_map["Heisen_Quasi_Sz_Sector_Ham"] = Heisen_Quasi_Sz_Sector_Ham_Para;
    model_para_read_map["Ising_Random_Simp_Cos_Ham"] = Ising_Random_Simp_Cos_Ham_Para;
    model_para_read_map["Ising_Quasi_Simp_Ham"] = Ising_Quasi_Simp_Ham_Para;
    model_para_read_map["Heisen_Modified_Random_Cos_Sz_Sector_Ham"] = Heisen_Modified_Random_Cos_Sz_Sector_Ham_Para;
    model_para_read_map["Heisen_Modified_Quasi_Sz_Sector_Ham"] = Heisen_Modified_Quasi_Sz_Sector_Ham_Para;
    model_para_read_map["Heisen_Con_Modified_Random_Cos_Sz_Sector_Ham"] = Heisen_Con_Modified_Random_Cos_Sz_Sector_Ham_Para;
    model_para_read_map["Heisen_Con_Modified_Quasi_Sz_Sector_Ham"] = Heisen_Con_Modified_Quasi_Sz_Sector_Ham_Para;

    task_para_read_map["Disorder_Transition"] = Disorder_Transition_Para;
    task_para_read_map["Single_Model_Time_Evolution"] = Single_Model_Time_Evolution_Para;
    task_para_read_map["Multi_Model_Time_Evolution"] = Multi_Model_Time_Evolution_Para;
    task_para_read_map["Op_Auto_Corr"] = Op_Auto_Corr_Para;
    task_para_read_map["Linked_Cluster"] = Linked_Cluster_Para;

    cout << "Read generic parameters:"<< endl;
    // Obtain generic parameters
    Generic_Para(parameters, count);

    cout << "Read output parameters:" << endl;
    // Obtain output parameters
    Output_Para(parameters, count);

    cout << "Read model parameters:" << endl;
    // Obtain model parameters
    map<string,para>::iterator it = model_para_read_map.find(parameters.generic.model);
    if (it != model_para_read_map.end()) it -> second(parameters, count);
    else{
        cout << "Data reading function for model " << parameters.generic.model << " cannot be found." << endl;
        abort();
    }

    cout << "Read task parameters:" << endl;
    // Obtain task parameters
    it = task_para_read_map.find(parameters.generic.task);
    if (it != task_para_read_map.end()) it -> second(parameters, count);
    else{
        cout << "Data reading function for task " << parameters.generic.task << " cannot be found." << endl;
        abort();
    }

    Eigen::initParallel();

    cout << "Start running." << endl;

    tasks_models.Task(parameters.generic.task)(parameters);

    cout << "Calculation finished." << endl;

    return 0;
}


