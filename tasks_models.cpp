//
// Created by Liangsheng Zhang on 4/15/15.
//

#include <iostream>
#include <utility>
#include <iostream>
#include "tasks_models.h"
#include "model_func.h"
#include "task_func.h"

using namespace std;

void TasksModels::Map_Construct_(){

    // Study transition property of a disorder model
    string task_name1;
    string task_type1;
    task_func task_function1;

    task_name1 = "Disorder_Transition";
    task_type1 = "All";
    task_function1 = &disorder_transition;
    Task_Map_Insert(task_name1, task_type1, task_function1);

    // Study time evolution of a single model
    string task_name2;
    string task_type2;
    task_func task_function2;

    task_name2 = "Single_Model_Time_Evolution";
    task_type2 = "All";
    task_function2 = &single_model_time_evolution;
    Task_Map_Insert(task_name2, task_type2, task_function2);

    // Study time evolution of a multiple identical models
    string task_name3;
    string task_type3;
    task_func task_function3;

    task_name3 = "Multi_Model_Time_Evolution";
    task_type3 = "All";
    task_function3 = &multi_model_time_evolution;
    Task_Map_Insert(task_name3, task_type3, task_function3);

    // Study the autocorrelation of some operator.
    // For now only under Floquet dynamics.
    string task_name4;
    string task_type4;
    task_func task_function4;

    task_name4 = "Op_Auto_Corr";
    task_type4 = "Floquet";
    task_function4 = &flo_op_auto_corr;
    Task_Map_Insert(task_name4, task_type4, task_function4);

    // Random Simple Ising Floquet Operator
    string model_name1;
    string model_type1;
    model_func model_function1;

    model_name1 = "Ising_Random_Simp_Flo";
    model_type1 = "Ising_Random_Simp_Floquet";
    model_function1 = &Flo_Evol_Ising_Random_Simp_Func;
    Model_Map_Insert(model_name1, model_type1, model_function1);

    // Random Simple Shift Real Ising Floquet Operator
    string model_name2;
    string model_type2;
    model_func model_function2;

    model_name2 = "Ising_Random_Simp_Shift_Real_Flo";
    model_type2 = "Ising_Random_Simp_Shift_Real_Floquet";
    model_function2 = &Flo_Evol_Ising_Random_Simp_Shift_Real_Func;
    Model_Map_Insert(model_name2, model_type2, model_function2);

    // Quasi-periodic Simple Shift Real Ising Floquet Operator
    string model_name3;
    string model_type3;
    model_func model_function3;

    model_name3 = "Ising_Quasi_Simp_Shift_Real_Flo";
    model_type3 = "Ising_Quasi_Simp_Shift_Real_Floquet";
    model_function3 = &Flo_Evol_Ising_Quasi_Simp_Shift_Real_Func;
    Model_Map_Insert(model_name3, model_type3, model_function3);

    // All Random Simple Shift Real Ising Floquet Operator
    string model_name4;
    string model_type4;
    model_func model_function4;

    model_name4 = "Ising_All_Random_Simp_Shift_Real_Flo";
    model_type4 = "Ising_All_Random_Simp_Shift_Real_Floquet";
    model_function4 = &Flo_Evol_Ising_All_Random_Simp_Shift_Real_Func;
    Model_Map_Insert(model_name4, model_type4, model_function4);

    // All Quasi-Periodic Simple Shift Real Ising Floquet Operator
    string model_name5;
    string model_type5;
    model_func model_function5;

    model_name5 = "Ising_All_Quasi_Simp_Shift_Real_Flo";
    model_type5 = "Ising_All_Quasi_Simp_Shift_Real_Floquet";
    model_function5 = &Flo_Evol_Ising_All_Quasi_Simp_Shift_Real_Func;
    Model_Map_Insert(model_name5, model_type5, model_function5);

    // Random Simple Shift Cos Real Ising Floquet Operator
    string model_name6;
    string model_type6;
    model_func model_function6;

    model_name6 = "Ising_Random_Simp_Shift_Cos_Real_Flo";
    model_type6 = "Ising_Random_Simp_Shift_Cos_Real_Floquet";
    model_function6 = &Flo_Evol_Ising_Random_Simp_Shift_Cos_Real_Func;
    Model_Map_Insert(model_name6, model_type6, model_function6);

    // Random Heisenberg Cos Total Z Spin Hamiltonian Operator
    string model_name7;
    string model_type7;
    model_func model_function7;

    model_name7 = "Heisen_Random_Cos_Sz_Sector_Ham";
    model_type7 = "Heisen_Random_Cos_Sz_Sector_Hamiltonian";
    model_function7 = &Ham_Evol_Heisen_Random_Cos_Sz_Sector_Func;
    Model_Map_Insert(model_name7, model_type7, model_function7);

    // Quasi-periodic Heisenberg Total Z Spin Hamiltonian Operator
    string model_name8;
    string model_type8;
    model_func model_function8;

    model_name8 = "Heisen_Quasi_Sz_Sector_Ham";
    model_type8 = "Heisen_Quasi_Sz_Sector_Hamiltonian";
    model_function8 = &Ham_Evol_Heisen_Quasi_Sz_Sector_Func;
    Model_Map_Insert(model_name8, model_type8, model_function8);

    // Random Simple Shift Cos Real Tau Ising Floquet Operator
    string model_name9;
    string model_type9;
    model_func model_function9;

    model_name9 = "Ising_Random_Simp_Shift_Cos_Real_Tau_Flo";
    model_type9 = "Ising_Random_Simp_Shift_Cos_Real_Tau_Floquet";
    model_function9 = &Flo_Evol_Ising_Random_Simp_Shift_Cos_Real_Tau_Func;
    Model_Map_Insert(model_name9, model_type9, model_function9);

    // Random Heisenberg Cos Sz Sector Shift Real Tau Floquet Operator
    string model_name10;
    string model_type10;
    model_func model_function10;

    model_name10 = "Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Flo";
    model_type10 = "Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Floquet";
    model_function10 = &Flo_Evol_Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Func;
    Model_Map_Insert(model_name10, model_type10, model_function10);

    // Quasi-periodic Heisenberg Sz Sector Shift Real Tau Floquet Operator
    string model_name11;
    string model_type11;
    model_func model_function11;

    model_name11 = "Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Flo";
    model_type11 = "Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Floquet";
    model_function11 = &Flo_Evol_Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Func;
    Model_Map_Insert(model_name11, model_type11, model_function11);
}

void TasksModels::Task_Map_Insert(const string& task_name, const string& task_type,
                                  const task_func& task_function){
    pair<string, task_func> task_type_op;

    if (tasks_.find(task_name) == tasks_.end()){ // No duplicate names
        task_type_op = make_pair(task_type, task_function);
        tasks_[task_name] = task_type_op;
    }
    else{
        cout << "Duplicate names for tasks." << endl;
        abort();
    }
}

void TasksModels::Model_Map_Insert(const string& model_name, const string& model_type,
                                   const model_func& model_function){
    pair<string, model_func> model_type_op;

    if (models_.find(model_name) == models_.end()){ // No duplicate names
        model_type_op = make_pair(model_type, model_function);
        models_[model_name] = model_type_op;
    }
    else{
        cout << "Duplicate names for models." << endl;
        abort();
    }
}

bool TasksModels::Task_Look_Up(const string& task_name) const {
    if (tasks_.find(task_name) == tasks_.end()){
        return false;
    }
    else return true;
}

bool TasksModels::Model_Look_Up(const string& model_name) const {
    if (models_.find(model_name) == models_.end()){
        return false;
    }
    else return true;
}

task_func TasksModels::Task(const string& task_name) {
    map<string, pair<string, task_func> >::const_iterator it;
    it = tasks_.find(task_name);
    if (it == tasks_.end()){
        cout << "The task desired is not found." << endl;
        cout << "Desired task: " << task_name << endl;
        Print_Task();
        abort();
    }
    else{
        task_type_ = it -> second.first;
        return it -> second.second;
    }
}

void TasksModels::Model(const string& model_name, const AllPara& parameters,
                        EvolOP*& model) const {
    map<string, pair<string, model_func> >::const_iterator it;
    it = models_.find(model_name);

    if (it == models_.end()){
        cout << "The model desired is not found." << endl;
        cout << "Desired model: " << model_name << endl;
        Print_Model();
        abort();
    }
    else{
        // Intrinsic type is constructed here
        string type;
        type = it -> second.second(parameters, model);

        if (it -> second.first != type){
            cout << "Model type is not consistent." << endl;
            cout << "Desired model type: " << it -> second.first << endl;
            cout << "Constructed model type: " << type << endl;
            abort();
        }

        // Assume the task_type_ belongs to the task which calls this model
        if (task_type_ != "All" && it -> second.first.find(task_type_) == string::npos){
            cout << "Model type and task type are not consistent." << endl;
            cout << "Model type: " << it -> second.first << endl;
            cout << "Task type: " << task_type_ << endl;
            abort();
        }
    }
}

void TasksModels::Print_Task() const {
    cout << "The tasks are: " << endl;
    map<string, pair<string, task_func> >::const_iterator it;

    for (it = tasks_.begin(); it != tasks_.end(); it++){
        cout << "Name: " << it -> first <<"  Type: " << it -> second.first << endl;
    }
    cout << endl;
}

void TasksModels::Print_Model() const {
    cout << "The models are: " << endl;
    map<string, pair<string, model_func> >::const_iterator it;

    for (it = models_.begin(); it != models_.end(); it++){
        cout << "Name: " << it -> first <<"  Type: " << it -> second.first << endl;
    }
    cout << endl;
}

