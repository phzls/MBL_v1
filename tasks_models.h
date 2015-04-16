//
// Created by Liangsheng Zhang on 4/15/15.
//

/*
 * This file contains a class which associates names with functions in task_func.h and
 * functions in model_func.h. Note it is assumed that all functions in task_func.h
 * just take parameter class as argument.
 */

#ifndef MBL_V1_TASK_MODEL_H
#define MBL_V1_TASK_MODEL_H

#include <map>
#include <utility>
#include <iostream>
#include <string>
#include "parameters.h"
#include "evol_op.h"

using namespace std;

typedef void (*task_func)(const AllPara&);
typedef string (*model_func)(const AllPara&, EvolOP*&);

class TasksModels
{
private:
    // Map for task names. The key is its name, and the value is the task type and
    // corresponding function call
    map<string, pair<string, task_func> > tasks_;

    // Map for model names. The key is its name, and the value is the model type, which
    // should be the same as the type obtained from the model class, and the corresponding
    // function call
    map<string, pair<string, model_func> > models_;

    // The recorded task type
    string task_type_;

    // Construct the two maps
    void Map_Construct_();

    // Insert task pair in map
    void Task_Map_Insert(const string&, const string&, const task_func&);

    // Insert model pair in map
    void Model_Map_Insert(const string&, const string&, const model_func&);

public:
    TasksModels() {Map_Construct_();}
    ~TasksModels();

    // Look up a name in tasks. If it exists, return true, otherwise false
    bool Task_Look_Up(const string&) const;

    // Look up a name in models. If it exists, return true, otherwise false
    bool Model_Look_Up(const string&) const;


    // According to the name of the task, return the task function pointer
    task_func Task(const string&);

    // According to the name of the model, return one model pointer
    void Model(const string&, const AllPara&, EvolOP*&) const;

    // Print all tasks names
    void Print_Task() const;

    // Print all models names
    void Print_Model() const;
};

#endif //MBL_V1_TASK_MODEL_H
