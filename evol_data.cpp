//
// Created by Liangsheng Zhang on 6/1/15.
//

#include <string>
#include <fstream>
#include <iomanip>
#include "evol_data.h"
#include "generic_func.h"

using namespace std;

/**
 ** This file contains some common functions in evol_data.h
 **/

EvolInfo::EvolInfo(const AllPara& parameters) {
    num_realization = parameters.generic.num_realizations;
    debug = parameters.generic.debug; // Whether print debug information
    threads_N = parameters.generic.threads_N; // Number of threads for parallelization
    time_step = parameters.evolution.time_step; // Number of time steps
    jump = parameters.evolution.jump; // jump of time points

    // Whether time changes logarithmically
    log_time = parameters.evolution.log_time;

    // The base under which time changes logarithmically
    log_time_jump = parameters.evolution.log_time_jump;

    left_size = parameters.evolution.left_size;
}

EvolData::EvolData(const AllPara& parameters): size_(parameters.generic.size), evol_info(parameters) {
    func_status_ = parameters.evolution.evol_compute;
    Data_Func_Map_Init_();

    if ( func_status_.size() != data_init_.size()){
        cout << "Number of registered functions in parameters for evolution is not the same as"
        << " number of registered functions that can be called." << endl;
        cout << "Functions in parameters:" << endl;
        for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
            cout << it -> first << endl;
        }
        cout << "Total Number: " << func_status_.size() << endl;

        cout << "Functions that can be called:" << endl;
        for (map<string, Data_Init>::iterator it = data_init_.begin(); it != data_init_.end();
             it ++){
            cout << it -> first << endl;
        }
        cout << "Total Number:" << data_init_.size() << endl;

        abort();
    }

    Name_Check_();

    // Initialize all quantities which will be computed
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second) ( this ->* (data_init_[it -> first]) ) (parameters);
    }
}

void EvolData::Print_All_Name() const{
    map<string, bool>::const_iterator it;

    for (it = func_status_.begin(); it != func_status_.end(); it ++){
        cout << it -> first << endl;
    }
}

void EvolData::Print_All_Status() const{
    map<string, bool>::const_iterator it;

    for (it = func_status_.begin(); it != func_status_.end(); it ++){
        cout << it -> first << "   ";

        if (it -> second) cout << "Would be computed." << endl;
        else cout << "Would not be computed." << endl;
    }
}

void EvolData::Name_Check_() const{
    map<string, bool>::const_iterator para_it;
    map<string, Data_Init>::const_iterator data_it;

    for (para_it = func_status_.begin(); para_it != func_status_.end(); para_it ++){
        data_it = data_init_.find(para_it -> first);
        if (data_it == data_init_.end()){
            cout << "Names in evolution are not consistent." << endl;

            cout << "Names in parameters:" << endl;
            map<string, bool>::const_iterator it;
            for (it = func_status_.begin(); it!= func_status_.end(); it++){
                cout << it -> first << endl;
            }

            cout << "Names in EvolData:" << endl;
            for (data_it = data_init_.begin(); data_it != data_init_.end(); data_it ++){
                cout << data_it -> first << endl;
            }
            abort();
        }
    }
}

void EvolData::Data_Compute(const VectorXcd& state, const StepInfo& info){
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second) ( this ->* (data_cal_[it -> first]) ) (state, info);
    }
}

void EvolData::Data_Compute(const MatrixXcd& state_density, const StepInfo& info){
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second) ( this ->* (data_cal_C_[it -> first]) ) (state_density, info);
    }
}

void EvolData::Init_Evol_Data(const InitEvolData& init_evol_data, const InitEvolInfo& init_evol_info) {
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second) ( this ->* (init_evol_[it -> first]) ) (init_evol_data, init_evol_info);
    }
}

void EvolData::Data_Output(const AllPara& parameters, const string& type_name){
    for (map<string, bool>::iterator it = func_status_.begin(); it != func_status_.end(); it++){
        if (it -> second){
            map<string, Data_Out>::iterator out_it = data_out_.find(it -> first);
            if (out_it != data_out_.end())
                ( this ->* (data_out_[it -> first]) ) (parameters, type_name);
        }
    }
}

void EvolData::Data_Func_Map_Init_(){
    map<string, Data_Init>::iterator init_it;
    map<string, Data_Cal>::iterator cal_it;
    map<string, Data_Cal_C>::iterator cal_C_it;
    map<string, Data_Out>::iterator out_it;
    map<string, Init_Evol>::iterator init_evol_it;

    // Entropy Per Model data
    string name1 = "Entropy_Per_Model";
    Data_Init init_func1 = &EvolData::Entropy_Per_Model_Init_;
    Init_Evol init_evol_func1 = &EvolData::Entropy_Per_Model_Evol_Init_;
    Data_Cal cal_func1 = &EvolData::Entropy_Per_Model_Cal_;
    Data_Cal_C cal_C_func1 = &EvolData::Entropy_Per_Model_Cal_C_;
    Data_Out out_func1 = &EvolData::Entropy_Per_Model_Out_;

    // Make sure the name has not been used before
    init_it = data_init_.find(name1);
    cal_it = data_cal_.find(name1);
    out_it = data_out_.find(name1);
    init_evol_it = init_evol_.find(name1);
    cal_C_it = data_cal_C_.find(name1);
    if (init_it != data_init_.end() || cal_it != data_cal_.end() || out_it != data_out_.end()
        || init_evol_it != init_evol_.end() || cal_C_it != data_cal_C_.end()){
        cout << name1 << " for evolution has appeared before." << endl;
        abort();
    }

    data_init_[name1] = init_func1;
    data_cal_[name1] = cal_func1;
    data_cal_C_[name1] = cal_C_func1;
    data_out_[name1] = out_func1;
    init_evol_[name1] = init_evol_func1;

    // Leftmost Spin Per Model data
    string name2 = "Leftmost_Spin_Per_Model";
    Data_Init init_func2 = &EvolData::Leftmost_Spin_Per_Model_Init_;
    Init_Evol init_evol_func2 = &EvolData::Leftmost_Spin_Per_Model_Evol_Init_;
    Data_Cal cal_func2 = &EvolData::Leftmost_Spin_Per_Model_Cal_;
    Data_Cal_C cal_C_func2 = &EvolData::Leftmost_Spin_Per_Model_Cal_C_;
    Data_Out out_func2 = &EvolData::Leftmost_Spin_Per_Model_Out_;

    // Make sure the name has not been used before
    init_it = data_init_.find(name2);
    cal_it = data_cal_.find(name2);
    out_it = data_out_.find(name2);
    init_evol_it = init_evol_.find(name2);
    cal_C_it = data_cal_C_.find(name2);
    if (init_it != data_init_.end() || cal_it != data_cal_.end() || out_it != data_out_.end()
        || init_evol_it != init_evol_.end() || cal_C_it != data_cal_C_.end()){
        cout << name2 << " for evolution has appeared before." << endl;
        abort();
    }

    data_init_[name2] = init_func2;
    data_cal_[name2] = cal_func2;
    data_cal_C_[name2] = cal_C_func2;
    data_out_[name2] = out_func2;
    init_evol_[name2] = init_evol_func2;

    // Check data_init_ and init_evol_ have the same size
    if (data_init_.size() != init_evol_.size()){
        cout << "Number of initializations in evolution is not the same as number of init_evol."
        << endl;
        cout << "Registered initializations:" << endl;
        for (init_it = data_init_.begin(); init_it != data_init_.end(); init_it ++){
            cout << init_it -> first << endl;
        }
        cout << "Total Number: " << data_init_.size();

        cout << "Registered init_evol:" << endl;
        for (init_evol_it = init_evol_.begin(); init_evol_it != init_evol_.end(); init_evol_it ++){
            cout << init_evol_it -> first << endl;
        }
        cout << "Total Number: " << init_evol_.size();
    }

    // Check data_init_ and data_cal_ have the same size
    if (data_init_.size() != data_cal_.size()){
        cout << "Number of initializations in evolution is not the same as number of calculations."
        << endl;
        cout << "Registered initializations:" << endl;
        for (init_it = data_init_.begin(); init_it != data_init_.end(); init_it ++){
            cout << init_it -> first << endl;
        }
        cout << "Total Number: " << data_init_.size();

        cout << "Registered calculations:" << endl;
        for (cal_it = data_cal_.begin(); cal_it != data_cal_.end(); cal_it ++){
            cout << cal_it -> first << endl;
        }
        cout << "Total Number: " << data_cal_.size();
    }

    // Check data_init_ and data_out_ have the same size
    if (data_init_.size() != data_out_.size() ){
        cout << "Number of initializations in evolution is not the same as number of output."
        << endl;
        cout << "Registered initializations:" << endl;
        for (init_it = data_init_.begin(); init_it != data_init_.end(); init_it ++){
            cout << init_it -> first << endl;
        }
        cout << "Total Number: " << data_init_.size();

        cout << "Registered output for each model:" << endl;
        for (out_it = data_out_.begin(); out_it != data_out_.end(); out_it ++){
            cout << out_it -> first << endl;
        }
    }
}

void EvolData::General_Output_(const AllPara& parameters, const vector<vector<vector<double> > >& data,
                               string filename) {
    const int time_step = parameters.evolution.time_step;
    const int width = parameters.output.width;
    const double step_size = parameters.evolution.step_size;
    const bool log_time = parameters.evolution.log_time;
    const int log_time_jump = parameters.evolution.log_time_jump;
    const bool markov_jump = parameters.evolution.markov_jump;
    const int markov_time_jump = parameters.evolution.markov_time_jump;
    const int jump = parameters.evolution.jump;

    ofstream fout( filename.c_str() );

    for (int t=0; t<time_step; t++){
        double time = t * step_size * jump; // An overflow still happens for entropy

        if (markov_jump) time *= markov_time_jump;

        if (log_time){
            long long int power = pow(log_time_jump,t);
            time = step_size * power;
        }

        if (parameters.evolution.sample_detail && (data.size() == 1)){
            fout << setw(10) << time;

            for (int i=0; i<data[0][t].size(); i++)
                fout << setw(width) << data[0][t][i];
            fout << endl;
        }
        else if (parameters.evolution.sample_detail){
            fout << setw(10) << time;
            for (int i=0; i < data.size(); i++){
                double mean, sd;
                generic_mean_sd(data[i][t], mean, sd);
                // For now errors are not outputted
                fout << setw(width) << mean;
            }
            fout << endl;
        }
        else{
            const int model_num = data.size();
            vector<double> mean(model_num);
            vector<double> sd(model_num);

            for (int i=0; i<model_num; i++){
                generic_mean_sd(data[i][t], mean[i], sd[i]);
            }

            double final_mean, final_sd;
            generic_mean_sd(mean, final_mean, final_sd);

            fout << setw(10) << time << setw(width) << final_mean << setw(width) << final_sd << endl;
        }
    }
}
