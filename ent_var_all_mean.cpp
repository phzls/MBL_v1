//
// Created by Liangsheng Zhang on 9/13/15.
//

#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <Eigen/Eigenvalues>
#include "disorder_model_transition.h"
#include "methods.h"
#include "generic_func.h"
#include "screen_output.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to ent_var_all_mean, which is the entropy "variance"
 ** for all eigenstates given a model. Here the mean used in calculation is the man from all
 ** realizations.
 ** The entropy is half chain entanglement entropy, so the spin chain must have even length.
 **/


/*
 * Initialize ent_var in model_data_ if it is not initialized in Ent_var_init_
 * The outer index is for J, and the inner index is for realization
 */
void DisorderModelTransition::Ent_var_all_mean_init_(AllPara const & parameters) {
    if (!flo_func_bool_map_["Entropy_Variance"]) Ent_var_init_(parameters);

    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.ent_mean.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.ent_mean[i].resize(num_realization);
        for (int j=0; j<num_realization;j++) model_data_.ent_mean[i][j] = 0;
    }

    ent_mean_keep_ = true;
}

/*
 * Compute ent_var if it is not computed in Ent_var_compute_
 */
void DisorderModelTransition::Ent_var_all_mean_compute_(AllPara const & parameters,
                                                        const EvolOP* model, const DisorderLocalInfo& local_info) {
    if (!flo_func_bool_map_["Entropy_Variance"]) {
        Ent_var_compute_(parameters, model, local_info);
        if (parameters.generic.debug) cout << "Ent_var_all_mean computes." << endl;
    }
    else if (parameters.generic.debug) cout << "Ent_var_all_mean does not compute." << endl;
}

/*
 * Output averages of mean entropy "variance" (over all eigenstates for each realization with single mean
 * for one J) and their standard deviations among all realizations for each J.
 * Note we only store "local" variance for each where the mean value is computed for each realization.
 * To get the desired result, we use the following formula:
 * Assume "local" mean for each realization is x, the mean of all x at one J for all realizations is y
 * "local" variance is w for each realization, then the "variance" with mean y for each realization is
 * x + N/(N-1)*(x-y)^z, where N is the number of eigenstates for each realization
 */
void DisorderModelTransition::Ent_var_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    const int dim = model_data_.model_dim; // Number of eigenstates for each realization

    if (model_data_.ent_var.size() != J_N){
        cout << "Not enough number of J for localentropy variance." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.ent_var.size() << endl;
        abort();
    }

    if (model_data_.ent_mean.size() != J_N){
        cout << "Not enough number of J for local entropy mean." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.ent_mean.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.ent_var[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for local entropy varaince." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.ent_var[i].size() << endl;
            abort();
        }
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.ent_mean[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for local entropy mean." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.ent_mean[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",entropy_variance_all_mean";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    vector<double> var(num_realizations); // Store re-calculated variace use single mean
    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;

        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        // Compute average across all realizations for one J, stored in mean
        generic_mean_sd(model_data_.ent_mean[i], mean, sd);

        // Compute entropy variance use single mean for each realization
        for(int j=0; j<num_realizations; j++){
            double local_mean = model_data_.ent_mean[i][j];
            var[i] = model_data_.ent_var[i][j] + double(dim)*(mean-local_mean)*(mean-local_mean)*(dim-1);
        }

        // Obtain average variance and its sd across all realizations for single J
        generic_mean_sd(var, mean, sd);

        fout << setw(10) << J << setw(width) << mean << setw(width) << sd << endl;
    }
}