//
// Created by Liangsheng Zhang on 9/14/15.
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
 ** This file contains implementation for Flo_level_stats_ave, average level statistics for
 ** floquet system.
 **
 ** The eval here are the phases of eigenvalues of floquet time evolution
 ** operator. At each level n, with ordered eval e_{n+1}, e_n and e_{n-1}, the level statistic
 ** considers s_n = e_{n+1} - e_n and s_{n-1} = e_n - e_{n-1}, and it is the minimum of the two
 ** divided by the maximum of the two.
 **
 ** Since eval are periodic, a periodic boundary condition is used, with proper 2*Pi shift of values
 ** across the boundary.
 **/

/*
 * Initialize flo_level_stat_ave in model_data_. The outer index is for J, and the inner index is for realization
 */
void DisorderModelTransition::Flo_level_stats_ave_init_(AllPara const & parameters) {
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
void DisorderModelTransition::Flo_level_stats_ave_compute_(AllPara const & parameters, const EvolOP* model,
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

    vector<double> phases(dim); // phases of all eigenvalues

    // Eigenvalues of Flouqet must all in general complex
    if (local_info.eval_type_real){
        cout << "Eigenvalues of Floquet operator cannot be real type." << endl;
        abort();
    }

    // the range of arg is between -Pi and Pi
    int index = 0;
    for(int i=0; i<local_info.eval_complex.size(); i++){
        for(int j=0; j<local_info.eval_complex[i].rows(); j++){
            if(index >= dim){
                cout << "Too many eigenvalues for dimension " << dim << endl;
                abort();
            }

            phases[index] = arg(local_info.eval_complex[i](j));
            index ++;
        }
    }
    if(index != dim-1){
        cout << "Incorrect number of eigenvalues" << endl;
        cout << "Expected: " << dim << endl;
        cout << "Obtained: " << index+1 << endl;
        abort();
    }

    // Sort phases in ascending order
    sort(phases.begin(), phases.end());

    double stats_ave = 0; // level statistics

    // Lower boundary case
    double s1 = phases[0] - (phases[dim-1] - 2*Pi);
    double s2 = phases[1] - phases[0];
    if(s1<s2) stats_ave = s1/s2;
    else stats_ave = s2/s1;

    for(int i=1; i<dim-1; i++){
        s1 = phases[i] - phases[i-1];
        s2 = phases[i+1] - phases[i];
        if(s1<s2) stats_ave += s1/s2;
        else stats_ave += s2/s1;
    }

    // Upper boundary case
    s1 = phases[dim-1] - phases[dim-2];
    s2 = phases[0] + 2*Pi - phases[dim-1];
    if(s1<s2) stats_ave += s1/s2;
    else stats_ave += s2/s1;

    stats_ave /= double(dim);

    model_data_.level_stat_ave[J][r] = stats_ave;

    if(debug){
        cout << "Realization: " << r << endl;
        cout << "Phases:" << endl;
        for(int i=0; i<phases.size(); i++) cout << phases[i] << endl;
        cout << "Average level statistic: " << stats_ave << endl;
    }
}

/*
 * Output averages of level statistics (which is an average over all eigenstates for each realization)
 * across all realizations for each J and the associated standard deviation
 */
void DisorderModelTransition::Flo_level_stats_ave_out_(AllPara const & parameters, const string& name) {
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
    << ",J_max=" << J_max << ",flo_level_statistic_ave";
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