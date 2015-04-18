//
// Created by Liangsheng Zhang on 4/17/15.
//

/*
 * Print out available tasks and models
 */

#include <iostream>
#include "tasks_models.h"
#include "randomc.h"

using namespace std;

CRandomMersenne RanGen_mersenne(time(NULL));

int EvolOP::model_num = 0;

TasksModels tasks_models;

int main(){

    tasks_models.Print_Model();

    tasks_models.Print_Task();

    return 0;
}


