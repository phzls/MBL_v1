//
// Created by Liangsheng Zhang on 9/24/15.
//

/*
 * This file implements functions in class OpCorrLocalPara
 */

#include "flo_op_auto_corr.h"

using namespace std;

OpCorrLocalPara::OpCorrLocalPara(const AllPara& parameters){
    // For now duplicate names are not checked
    string name = "tau_choice";
    para_update_func_[name] = &OpCorrLocalPara::Tau_update_;
    para_init_func_[name] = &OpCorrLocalPara::Tau_init_;
    para_update_bool_[name] = parameters.flo_op_auto_corr.tau_choice;

    // Check that it has exactly one true
    int true_count = 0;
    map<string, bool>::iterator it;
    for(it = para_update_bool_.begin(); it != para_update_bool_.end(); it++){
        if(it -> second){
            ++ true_count;
            para_name_ = it -> first;
        }
    }

    if(true_count != 1){
        cout << "The number of parameters to update in flo_op_auto_corr is not correct." << endl;
        cout << "1 is expected. Parameters required to update are:" << endl;
        for(it = para_update_bool_.begin(); it != para_update_bool_.end(); it ++){
            if(it -> second) cout << it -> first << endl;
        }
        abort();
    }

    if(para_update_func_.size() != para_update_bool_.size()){
        cout << "Incompatible size for bool and updating functions in flo_corr_local_para" << endl;
        cout << "In bool:" << endl;
        for(it = para_update_bool_.begin(); it != para_update_bool_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << para_update_bool_.size() << endl;
        cout << "In func:" << endl;
        map<string, Para_cal>::iterator it2;
        for(it2 = para_update_func_.begin(); it2 != para_update_func_.end(); it2 ++){
            cout << it2 -> first << endl;
        }
        cout << "Total Number: " << para_update_func_.size() << endl;
    }

    if(para_init_func_.size() != para_update_bool_.size()){
        cout << "Incompatible size for bool and initialization functions in flo_corr_local_para" << endl;
        cout << "In bool:" << endl;
        for(it = para_update_bool_.begin(); it != para_update_bool_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << para_update_bool_.size() << endl;
        cout << "In func:" << endl;
        map<string, Para_init>::iterator it2;
        for(it2 = para_init_func_.begin(); it2 != para_init_func_.end(); it2 ++){
            cout << it2 -> first << endl;
        }
        cout << "Total Number: " << para_init_func_.size() << endl;
    }

    // Update parameters
    ( this ->*(para_init_func_[para_name_]) )(parameters);

    para_pts_ = parameters.flo_op_auto_corr.para_pts;
    if(para_pts_ < 1){
        cout << "At least one instance of parameters must be allowed." << endl;
        cout << "Obtained number of instance: " << para_pts_ << endl;
        abort();
    }
}

void OpCorrLocalPara::Tau_init_(const AllPara& parameters) {
    min_para_ = parameters.flo_op_auto_corr.tau_min;
    max_para_ = parameters.flo_op_auto_corr.tau_max;
    para_ = "tau";
}

double OpCorrLocalPara::Tau_update_(int n, AllPara& local_para) {
    if(para_pts_ == 1){
        if(n>0){
            cout << "There is only one instance of parameters." << endl;
            cout << "Current step: " << n << endl;
            abort();
        }
        local_para.floquet.tau = min_para_;
        return local_para.floquet.tau;
    }

    double space = (max_para_ - min_para_) / double(para_pts_ - 1);
    local_para.floquet.tau = min_para_ + n * space;
    return local_para.floquet.tau;
}

