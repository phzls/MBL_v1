//
// Created by Liangsheng Zhang on 4/15/15.
//

/*
 * This file contains functions related to zz_corr, which is the end-to-end zz correlation square used in studying
 * floquet model transition from thermalization to localization
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
 * Initialize zz_corr_square in model_data_. The outer index is for J, and the inner index is for realization
 */
void DisorderModelTransition::ZZ_corr_square_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.zz_corr_square.resize(J_N);

    for (int i=0; i<J_N; i++){
    model_data_.zz_corr_square[i].resize(num_realization);

    for (int j=0; j<num_realization;j++) model_data_.zz_corr_square[i][j] = 0;
    }

    op_diag_ = true;
    op_dial_evec_keep_ = true;
}

/*
 * Compute zz_corr_square given J index and realization index
 */
void DisorderModelTransition::ZZ_corr_square_compute_(AllPara const & parameters,
                                                      const EvolOP* model, const DisorderLocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    const int dim = model -> Get_Dim();
    const int half_dim = dim / 2; // Below which leftmost spin is -1

    // Record values for all eigenstate
    double zz_square = 0;

    if (local_info.evec_type_real){
        // Vector for eigenvectors in basic binary basis
        vector<vector<double > > evec_basic(dim);
        for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, local_info.evec_real, evec_basic);

        for (int i=0; i<evec_basic.size(); i++){
            // Record values for each eigenstate
            double left_temp = 0;
            double right_temp = 0;
            double left_right_temp = 0;

            for (int j=0; j<evec_basic[i].size(); j++){
                int right_spin = -1;
                int left_spin = 1;
                if (j & 1 == 1) right_spin = 1;
                if (j < half_dim ) left_spin = -1;

                left_temp += left_spin * evec_basic[i][j] * evec_basic[i][j];
                right_temp += right_spin * evec_basic[i][j] * evec_basic[i][j];
                left_right_temp += left_spin * right_spin * evec_basic[i][j] * evec_basic[i][j];
            }

            zz_square += (left_right_temp - left_temp*right_temp ) * (left_right_temp - left_temp*right_temp );
            if (debug){
                cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
                cout << "Average left spin: " << left_temp << endl;
                cout << "Average right spin: " << right_temp << endl;
                cout << "Average left X right: " << left_right_temp << endl;
            }
        }
    }
    else{
        // Vector for eigenvectors in basic binary basis
        vector<vector<complex<double> > > evec_basic(dim);
        for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, local_info.evec_complex, evec_basic);

        for (int i=0; i<evec_basic.size(); i++){
            // Record values for each eigenstate
            double left_temp = 0;
            double right_temp = 0;
            double left_right_temp = 0;

            for (int j=0; j<evec_basic[i].size(); j++){
                int right_spin = -1;
                int left_spin = 1;
                if (j & 1 == 1) right_spin = 1;
                if (j < half_dim ) left_spin = -1;

                left_temp += left_spin * norm(evec_basic[i][j]);
                right_temp += right_spin * norm(evec_basic[i][j]);
                left_right_temp += left_spin * right_spin * norm(evec_basic[i][j]);
            }

            zz_square += (left_right_temp - left_temp*right_temp ) * (left_right_temp - left_temp*right_temp );
            if (debug){
                cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
                cout << "Average left spin: " << left_temp << endl;
                cout << "Average right spin: " << right_temp << endl;
                cout << "Average left X right: " << left_right_temp << endl;
            }
        }
    }

    model_data_.zz_corr_square[local_info.J_index][local_info.realization_index] += zz_square / double(dim);
    if (debug){
        cout << "Realization " << local_info.realization_index << " average correlation square: "
        << zz_square / double(dim) << endl;
        cout << endl;
    }
}

/*
 * Output averages of mean end-to-end zz correlation squares (over all eigenstates for each realization) and their
 * standard deviations among all realizations for each J
 */
void DisorderModelTransition::ZZ_corr_square_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.zz_corr_square.size() != J_N){
        cout << "Not enough number of J for zz correlation square." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.zz_corr_square.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.zz_corr_square[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for zz correlation square." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.zz_corr_square[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
             << ",J_max=" << J_max << ",zz_corr_square";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        generic_mean_sd(model_data_.zz_corr_square[i], mean, sd);
        fout << setw(10) << J << setw(width) << mean << setw(width) << sd << endl;
    }
}


