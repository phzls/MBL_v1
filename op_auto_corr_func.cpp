//
// Created by Liangsheng Zhang on 9/24/15.
//

/*
 * This file implements the class OpAutoCorr
 */

#include "flo_op_auto_corr.h"

using namespace std;

OpAutoCorr::OpAutoCorr(const AllPara& parameters) {
    map<string, Op_init >::iterator init_it;
    map<string, Op_set_up>::iterator set_up_it;
    map<string, Op_corr_cal >::iterator cal_it;
    map<string, Op_corr_out >::iterator out_it;

    op_auto_corr_map_.insert(parameters.flo_op_auto_corr.op_corr_map.begin(),
                              parameters.flo_op_auto_corr.op_corr_map.end());

    // Energy autocorrelation
    string name1 = "Energy";
    Op_init init_func1 = &OpAutoCorr::Energy_corr_init_;
    Op_set_up set_up_func1 = &OpAutoCorr::Energy_corr_set_up_;
    Op_corr_cal cal_func1 = &OpAutoCorr::Energy_corr_compute_;
    Op_corr_out out_func1 = &OpAutoCorr::Energy_corr_out_;

    // Make sure the name has not been used before
    init_it = op_init_map_.find(name1);
    set_up_it = op_set_up_map_.find(name1);
    cal_it = op_corr_cal_map_.find(name1);
    out_it = op_corr_out_map_.find(name1);
    if (init_it != op_init_map_.end() || set_up_it != op_set_up_map_.end() ||
        cal_it != op_corr_cal_map_.end() || out_it != op_corr_out_map_.end()){
        cout << name1 << " for op_auto_corr has appeared before." << endl;
        abort();
    }

    op_init_map_[name1] = init_func1;
    op_set_up_map_[name1] = set_up_func1;
    op_corr_cal_map_[name1] = cal_func1;
    op_corr_out_map_[name1] = out_func1;

    // Check the number of function
    if ( op_init_map_.size() != op_corr_cal_map_.size() ){

        cout << "Number of initializing functions in op_auto_corr is not the same as number of registered functions"
        << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, Op_init >::iterator it = op_init_map_.begin();
             it != op_init_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << op_init_map_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Op_corr_cal >::iterator it = op_corr_cal_map_.begin();
             it != op_corr_cal_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << op_corr_cal_map_.size() << endl;

        abort();
    }

    if ( op_set_up_map_.size() != op_corr_cal_map_.size() ){

        cout << "Number of setting up functions in op_auto_corr is not the same as number of registered functions"
        << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, Op_set_up >::iterator it = op_set_up_map_.begin();
             it != op_set_up_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << op_set_up_map_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Op_corr_cal >::iterator it = op_corr_cal_map_.begin();
             it != op_corr_cal_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << op_corr_cal_map_.size() << endl;

        abort();
    }

    if ( op_corr_out_map_.size() != op_corr_cal_map_.size() ){

        cout << "Number of output functions in op_auto_corr is not the same as number of registered functions"
        << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, Op_corr_out>::iterator it = op_corr_out_map_.begin();
             it != op_corr_out_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << op_corr_out_map_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Op_corr_cal >::iterator it = op_corr_cal_map_.begin();
             it != op_corr_cal_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << op_corr_cal_map_.size() << endl;

        abort();
    }

    if ( op_auto_corr_map_.size() != op_corr_cal_map_.size() ){

        cout << "Number of registered functions in op_auto_corr is different from the number of "
        << "registered functions that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, bool>::iterator it = op_auto_corr_map_.begin();
             it != op_auto_corr_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << op_auto_corr_map_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Op_corr_cal>::iterator it = op_corr_cal_map_.begin();
             it != op_corr_cal_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << op_corr_cal_map_.size() << endl;

        abort();
    }

    // Initialize data structures
    for (map<string, bool>::iterator it = op_auto_corr_map_.begin(); it != op_auto_corr_map_.end(); it++){
        if (it -> second) ( this ->*(op_init_map_[it -> first]) ) (parameters);
    }
}

void OpAutoCorr::SetUp(const AllPara& parameters, const EvolOP* evol) {
    // Set up the operators used
    for (map<string, bool>::iterator it = op_auto_corr_map_.begin(); it != op_auto_corr_map_.end(); it++){
        if (it -> second) ( this ->*(op_set_up_map_[it -> first]) ) (parameters, evol);
    }
}

void OpAutoCorr::Compute(AllPara const & parameters, OpCorrLocalInfo const & local_info) {
    // Calculation
    for (map<string, bool>::iterator it = op_auto_corr_map_.begin(); it != op_auto_corr_map_.end(); it++){
        if (it -> second) ( this ->*(op_corr_cal_map_[it -> first]) )(parameters, local_info);
    }
}

void OpAutoCorr::Output(const AllPara& parameters, const OpCorrLocalPara& op_corr_local_para, const string& name) {
    // Output file
    for (map<string, bool>::iterator it = op_auto_corr_map_.begin(); it != op_auto_corr_map_.end(); it++){
        if (it -> second) ( this ->*(op_corr_out_map_[it -> first]) )(parameters, op_corr_local_para, name);
    }
}
