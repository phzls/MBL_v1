//
// Created by Liangsheng Zhang on 4/15/15.
//

#include <iostream>
#include <cstdlib>
#include "disorder_model_transition.h"

using namespace std;

void DisorderModelTransition::map_initialize_(const AllPara& parameters) {

    map<string, Flo_init >::iterator init_it;
    map<string, Flo_func >::iterator cal_it;
    map<string, Flo_out >::iterator out_it;

    flo_func_bool_map_.insert(parameters.transition.flo_transition_compute.begin(),
    parameters.transition.flo_transition_compute.end());
    /*for (map<string,bool>::const_iterator para_it = parameters.transition.flo_transition_compute.begin();
            para_it != parameters.transition.flo_transition_compute.end(); para_it++){

    }*/

    // end-to-end z-z correlation square
    string name1 = "ZZ_Correlation_Square";
    Flo_init init_func1 = &DisorderModelTransition::ZZ_corr_square_init_;
    Flo_func cal_func1 = &DisorderModelTransition::ZZ_corr_square_compute_;
    Flo_out out_func1 = &DisorderModelTransition::ZZ_corr_square_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name1);
    cal_it = flo_func_map_.find(name1);
    out_it = flo_out_map_.find(name1);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name1 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name1] = init_func1;
    flo_func_map_[name1] = cal_func1;
    flo_out_map_[name1] = out_func1;

    // entropy variance
    string name2 = "Entropy_Variance";
    Flo_init init_func2 = &DisorderModelTransition::Ent_var_init_;
    Flo_func cal_func2 = &DisorderModelTransition::Ent_var_compute_;
    Flo_out out_func2 = &DisorderModelTransition::Ent_var_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name2);
    cal_it = flo_func_map_.find(name2);
    out_it = flo_out_map_.find(name2);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name2 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name2] = init_func2;
    flo_func_map_[name2] = cal_func2;
    flo_out_map_[name2] = out_func2;

    // end-to-end z-z time correlation square
    string name3 = "ZZ_Time_Correlation";
    Flo_init init_func3 = &DisorderModelTransition::ZZ_time_corr_init_;
    Flo_func cal_func3 = &DisorderModelTransition::ZZ_time_corr_compute_;
    Flo_out out_func3 = &DisorderModelTransition::ZZ_time_corr_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name3);
    cal_it = flo_func_map_.find(name3);
    out_it = flo_out_map_.find(name3);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name3 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name3] = init_func3;
    flo_func_map_[name3] = cal_func3;
    flo_out_map_[name3] = out_func3;

    // end-to-end z-z time correlation components
    string name4 = "ZZ_Time_Correlation_Components";
    Flo_init init_func4 = &DisorderModelTransition::ZZ_time_corr_component_init_;
    Flo_func cal_func4 = &DisorderModelTransition::ZZ_time_corr_component_compute_;
    Flo_out out_func4 = &DisorderModelTransition::ZZ_time_corr_component_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name4);
    cal_it = flo_func_map_.find(name4);
    out_it = flo_out_map_.find(name4);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name4 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name4] = init_func4;
    flo_func_map_[name4] = cal_func4;
    flo_out_map_[name4] = out_func4;

    // z-z correlation square at all distances
    string name5 = "ZZ_All_Correlation_Square";
    Flo_init init_func5 = &DisorderModelTransition::ZZ_all_corr_square_init_;
    Flo_func cal_func5 = &DisorderModelTransition::ZZ_all_corr_square_compute_;
    Flo_out out_func5 = &DisorderModelTransition::ZZ_all_corr_square_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name5);
    cal_it = flo_func_map_.find(name5);
    out_it = flo_out_map_.find(name5);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name5 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name5] = init_func5;
    flo_func_map_[name5] = cal_func5;
    flo_out_map_[name5] = out_func5;

    // end-to-end z-z correlation square for all samples
    string name6 = "ZZ_Correlation_Square_All_Sample";
    Flo_init init_func6 = &DisorderModelTransition::ZZ_corr_square_all_sample_init_;
    Flo_func cal_func6 = &DisorderModelTransition::ZZ_corr_square_all_sample_compute_;
    Flo_out out_func6 = &DisorderModelTransition::ZZ_corr_square_all_sample_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name6);
    cal_it = flo_func_map_.find(name6);
    out_it = flo_out_map_.find(name6);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name6 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name6] = init_func6;
    flo_func_map_[name6] = cal_func6;
    flo_out_map_[name6] = out_func6;

    // end-to-end z-z time correlation all sample
    string name7 = "ZZ_Time_Correlation_All_Sample";
    Flo_init init_func7 = &DisorderModelTransition::ZZ_time_corr_all_sample_init_;
    Flo_func cal_func7 = &DisorderModelTransition::ZZ_time_corr_all_sample_compute_;
    Flo_out out_func7 = &DisorderModelTransition::ZZ_time_corr_all_sample_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name7);
    cal_it = flo_func_map_.find(name7);
    out_it = flo_out_map_.find(name7);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name7 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name7] = init_func7;
    flo_func_map_[name7] = cal_func7;
    flo_out_map_[name7] = out_func7;

    // log z-z correlation square
    string name8 = "Log_ZZ_Correlation_Square";
    Flo_init init_func8 = &DisorderModelTransition::Log_ZZ_corr_square_init_;
    Flo_func cal_func8 = &DisorderModelTransition::Log_ZZ_corr_square_compute_;
    Flo_out out_func8 = &DisorderModelTransition::Log_ZZ_corr_square_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name8);
    cal_it = flo_func_map_.find(name8);
    out_it = flo_out_map_.find(name8);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name8 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name8] = init_func8;
    flo_func_map_[name8] = cal_func8;
    flo_out_map_[name8] = out_func8;

    // log z-z correlation square at all distances
    string name9 = "Log_ZZ_All_Correlation_Square";
    Flo_init init_func9 = &DisorderModelTransition::Log_ZZ_all_corr_square_init_;
    Flo_func cal_func9 = &DisorderModelTransition::Log_ZZ_all_corr_square_compute_;
    Flo_out out_func9 = &DisorderModelTransition::Log_ZZ_all_corr_square_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name9);
    cal_it = flo_func_map_.find(name9);
    out_it = flo_out_map_.find(name9);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name9 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name9] = init_func9;
    flo_func_map_[name9] = cal_func9;
    flo_out_map_[name9] = out_func9;

    // z-z time correlation at all distances
    string name10 = "ZZ_All_Time_Correlation";
    Flo_init init_func10 = &DisorderModelTransition::ZZ_all_time_corr_init_;
    Flo_func cal_func10 = &DisorderModelTransition::ZZ_all_time_corr_compute_;
    Flo_out out_func10 = &DisorderModelTransition::ZZ_all_time_corr_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name10);
    cal_it = flo_func_map_.find(name10);
    out_it = flo_out_map_.find(name10);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name10 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name10] = init_func10;
    flo_func_map_[name10] = cal_func10;
    flo_out_map_[name10] = out_func10;

    // entropy "variance" using one mean at each J
    string name11 = "Entropy_Variance_All_Mean";
    Flo_init init_func11 = &DisorderModelTransition::Ent_var_all_mean_init_;
    Flo_func cal_func11 = &DisorderModelTransition::Ent_var_all_mean_compute_;
    Flo_out out_func11 = &DisorderModelTransition::Ent_var_all_mean_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name11);
    cal_it = flo_func_map_.find(name11);
    out_it = flo_out_map_.find(name11);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name11 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name11] = init_func11;
    flo_func_map_[name11] = cal_func11;
    flo_out_map_[name11] = out_func11;

    // Average level statistics for Floquet system
    string name12 = "Flo_Level_Stats_Ave";
    Flo_init init_func12 = &DisorderModelTransition::Flo_level_stats_ave_init_;
    Flo_func cal_func12 = &DisorderModelTransition::Flo_level_stats_ave_compute_;
    Flo_out out_func12 = &DisorderModelTransition::Flo_level_stats_ave_out_;

    // Make sure the name has not been used before
    init_it = flo_init_map_.find(name12);
    cal_it = flo_func_map_.find(name12);
    out_it = flo_out_map_.find(name12);
    if (init_it != flo_init_map_.end() || cal_it != flo_func_map_.end() || out_it != flo_out_map_.end()){
        cout << name12 << " for disorder transition has appeared before." << endl;
        abort();
    }

    flo_init_map_[name12] = init_func12;
    flo_func_map_[name12] = cal_func12;
    flo_out_map_[name12] = out_func12;

    // Check the number of function
    if ( flo_init_map_.size() != flo_func_map_.size() ){

        cout << "Number of initializing functions in flo_transition is not the same as number of registered functions"
        << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, Flo_init >::iterator it = flo_init_map_.begin();
             it != flo_init_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << flo_init_map_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Flo_func>::iterator it = flo_func_map_.begin();
             it != flo_func_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << flo_func_map_.size() << endl;

        abort();
    }

    if ( flo_out_map_.size() != flo_func_map_.size() ){

        cout << "Number of output functions in flo_transition is not the same as number of registered functions"
        << "that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, Flo_out>::iterator it = flo_out_map_.begin();
             it != flo_out_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << flo_out_map_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Flo_func>::iterator it = flo_func_map_.begin();
             it != flo_func_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << flo_func_map_.size() << endl;

        abort();
    }

    if ( flo_func_bool_map_.size() != flo_func_map_.size() ){

        cout << "Number of registered functions in flo_transition is different from the number of "
        << "registered functions that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, bool>::iterator it = flo_func_bool_map_.begin();
             it != flo_func_bool_map_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << flo_func_bool_map_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Flo_func>::iterator it = flo_func_map_.begin();
             it != flo_func_map_.end(); it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << flo_func_map_.size() << endl;

        abort();
    }

    // Initialize data structures
    for (map<string, bool>::iterator it = flo_func_bool_map_.begin(); it != flo_func_bool_map_.end(); it++){
        if (it -> second) ( this ->*(flo_init_map_[it -> first]) ) (parameters);
    }
}

void DisorderModelTransition::Compute(AllPara const & parameters, EvolOP const * floquet,
                                      DisorderLocalInfo const & local_info) {

    for (map<string, bool>::iterator it = flo_func_bool_map_.begin(); it != flo_func_bool_map_.end(); it++){
        if (it -> second) ( this ->*(flo_func_map_[it -> first]) )(parameters, floquet, local_info);
    }
}

void DisorderModelTransition::Output(AllPara const & parameters, const string& name) {
    for (map<string, bool>::iterator it = flo_func_bool_map_.begin(); it != flo_func_bool_map_.end(); it++){
        if (it -> second) ( this ->*(flo_out_map_[it -> first]) )(parameters, name);
    }
}


