//
// Created by Liangsheng Zhang on 4/15/15.
//

/*
 * This file includes all the tasks that can be called from the main function.
 */

#ifndef MBL_V1_TASK_FUNC_H
#define MBL_V1_TASK_FUNC_H

#include "parameters.h"

/*
 * Study transition properties of a disorder model
 */
void disorder_transition(const AllPara&);

/*
 * Study time evolution under a single model
 */
void single_model_time_evolution(const AllPara&);

/*
 * Study time evolution with multiple identical models
 */
void multi_model_time_evolution(const AllPara&);

/*
 * Study the autocorrelation of some operator.
 * For now only under Floquet dynamics.
 */
void flo_op_auto_corr(const AllPara&);

/*
 * Study the coefficients of linked clusters.
 * For now all model parameters must be real
 */
void linked_cluster(const AllPara&);

#endif //MBL_V1_TASK_FUNC_H
