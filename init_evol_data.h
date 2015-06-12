//
// Created by Liangsheng Zhang on 6/12/15.
//

#ifndef MBL_V1_INIT_EVOL_DATA_H
#define MBL_V1_INIT_EVOL_DATA_H

/*
 * This file contains data from initial condition which will be used later for time evolution
 */
struct InitEvolData{
    double infinite_time_leftmost_spin;
};

struct InitEvolInfo{
    int model;
    int realization;
    bool debug;
};

#endif //MBL_V1_INIT_EVOL_DATA_H
