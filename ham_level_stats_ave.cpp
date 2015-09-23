//
// Created by Liangsheng Zhang on 9/21/15.
//

#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <Eigen/Eigenvalues>
#include "disorder_model_transition.h"
#include "methods.h"
#include "generic_func.h"
#include "screen_output.h"
#include "constants.h"

/**
 ** This file contains implementation for Ham_level_stats_ave, average level statistics for
 ** floquet system.
 **
 ** The eval here are the eigenvalues of Hamiltonian operator
 ** At each level n, with ordered eval e_{n+1}, e_n and e_{n-1}, the level statistic
 ** considers s_n = e_{n+1} - e_n and s_{n-1} = e_n - e_{n-1}, and it is the minimum of the two
 ** divided by the maximum of the two.
 **
 ** Since we are working with the mid_half of the spectrum, the one above and the one below
 ** are used for boundaries
 **/

/*
 * Initialize ham_level_stat_ave in model_data_. The outer index is for J, and the inner index is for realization
 */
void DisorderModelTransition::Ham_level_stats_ave_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.level_stat_ave.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.level_stat_ave[i].resize(num_realization);

        for (int j=0; j<num_realization;j++) model_data_.level_stat_ave[i][j] = 0;
    }

    op_diag_ = true;
}

/*
 * Compute levels statistics given J index and realization index
 */
void DisorderModelTransition::Ham_level_stats_ave_compute_(AllPara const & parameters, const EvolOP* model,
                                                           const DisorderLocalInfo& local_info) {

    const bool debug = parameters.generic.debug;
    const int dim = model -> Get_Dim();

    const int J = local_info.J_index;
    const int r = local_info.realization_index;

    if (dim < 3){
        cout << "There must be at least 3 eigenvalues" << endl;
        cout << "System dimension: " << dim << endl;
        abort();
    }

    // Eigenvalues of Hamiltonian must all be real
    if (!local_info.eval_type_real){
        cout << "Eigenvalues of Hamiltonian must be real type." << endl;
        abort();
    }

    // For now the program only works with 1 sector, as it is not
    // clear how to pool eigenvalues of different sectors together
    if(local_info.eval_real.size() != 1){
        cout << "ham_level_stats_ave only works with 1 sector for now." << endl;
        cout << "Obtained number of sectors: " << local_info.eval_real.size() << endl;
        abort();
    }

    double stats_ave = 0; // level statistics

    if(local_info.mid_half_spectrum){ // Only middle half of the spectrum is used
        int mid_half_size = local_info.eval_real[0].size();

        if(mid_half_size != dim-2*(dim/4)){
            cout << "Incorrect number of middle half eigenvalues" << endl;
            cout << "Expected: " << dim - 2*(dim/4) << endl;
            cout << "Obtained: " << mid_half_size << endl;
            abort();
        }

        // Lower boundary case
        double s1 = local_info.eval_real[0][0] - local_info.lower_eval;
        double s2 = local_info.upper_eval - local_info.eval_real[0][0];
        if(mid_half_size > 1) s2 = local_info.eval_real[0][1] - local_info.eval_real[0][0];
        if(s1<s2) stats_ave = s1/s2;
        else stats_ave = s2/s1;

        for(int i=1; i<mid_half_size-1; i++){
            s1 = local_info.eval_real[0][i] - local_info.eval_real[0][i-1];
            s2 = local_info.eval_real[0][i+1] - local_info.eval_real[0][i];
            if(s1<s2) stats_ave += s1/s2;
            else stats_ave += s2/s1;
        }

        // Upper boundary case
        if(mid_half_size>1){
            s1 = local_info.eval_real[0][mid_half_size-1] - local_info.eval_real[0][mid_half_size-2];
            s2 = local_info.upper_eval - local_info.eval_real[0][mid_half_size-1];
            if(s1<s2) stats_ave += s1/s2;
            else stats_ave += s2/s1;
        }

        stats_ave /= double(mid_half_size);
    }
    else{ // All eigenvalues are taken
        if( local_info.eval_real[0].size() != dim ){
            cout << "Incorrect number of eigenvalues for ham_level_stats" << endl;
            cout << "Expected: " << dim << endl;
            cout << "Obtained: " << local_info.eval_real[0].size() << endl;
            abort();
        }

        for(int i=1; i<dim-1; i++){
            double s1 = local_info.eval_real[0][i] - local_info.eval_real[0][i-1];
            double s2 = local_info.eval_real[0][i+1] - local_info.eval_real[0][i];
            if(s1<s2) stats_ave += s1/s2;
            else stats_ave += s2/s1;
        }

        stats_ave /= double(dim-2);
    }

    model_data_.level_stat_ave[J][r] = stats_ave;

    if(debug){
        cout << "Realization: " << r << endl;
        cout << "Average level statistic: " << stats_ave << endl;
    }
}

/*
 * Output averages of level statistics (which is an average over all eigenstates for each realization)
 * across all realizations for each J and the associated standard deviation
 */
void DisorderModelTransition::Ham_level_stats_ave_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.level_stat_ave.size() != J_N){
        cout << "Not enough number of J for entropy variance." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.level_stat_ave.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.level_stat_ave[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for entropy varaince." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.level_stat_ave[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",ham_level_statistic_ave";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        generic_mean_sd(model_data_.level_stat_ave[i], mean, sd);
        fout << setw(10) << J << setw(width) << mean << setw(width) << sd << endl;
    }
}

