//
// Created by Liangsheng Zhang on 4/23/15.
//

/*
 * This file contains functions which output zz_time_corr, end-to-end zz time four point correlation
 * used in studying floquet model transition from thermalization to localization, from all realizations.
 */

#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include "disorder_model_transition.h"
#include "methods.h"
#include "generic_func.h"

using namespace std;

/*
 * Initialize zz_time_corr in model_data_ if it is not initialized in ZZ_time_corr_all_sample_init_
 */
void DisorderModelTransition::ZZ_time_corr_all_sample_init_(AllPara const & parameters) {
    if (!flo_func_bool_map_["ZZ_Time_Correlation"]) {
        ZZ_time_corr_init_(parameters);
        if (parameters.generic.debug) cout << "ZZ_time_corr_all_sample initializes." << endl;
    }
    else if (parameters.generic.debug) cout << "ZZ_time_corr_all_sample does not initialize." << endl;
}

/*
 * Compute zz_time_corr if it is not computed in ZZ_time_corr_compute_
 */
void DisorderModelTransition::ZZ_time_corr_all_sample_compute_(AllPara const & parameters,
                                                    const EvolOP* model,
                                                    const DisorderLocalInfo& local_info) {

    if (!flo_func_bool_map_["ZZ_Time_Correlation"]) {
        ZZ_time_corr_compute_(parameters, model, local_info);
        if (parameters.generic.debug) cout << "ZZ_time_corr_all_sample computes." << endl;
    }
    else if (parameters.generic.debug) cout << "ZZ_time_corr_all_sample does not compute." << endl;
}

/*
 * Output averages of mean end-to-end zz time four-point correlation squares and their standard deviations for all
 * realizations for each J
 */
void DisorderModelTransition::ZZ_time_corr_all_sample_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.zz_time_corr.size() != J_N){
        cout << "Not enough number of J for zz time correlation." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.zz_time_corr.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.zz_time_corr[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for zz time correlation." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.zz_time_corr[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",zz_time_corr_all_sample";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        fout << setw(10) << J;

        for (int j=0; j<model_data_.zz_time_corr[i].size();j++){
            fout << setw(width) << model_data_.zz_time_corr[i][j];
        }
        fout << endl;
    }
}





