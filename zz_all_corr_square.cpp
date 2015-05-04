//
// Created by Liangsheng Zhang on 4/15/15.
//

/*
 * This file contains functions related to zz_all_corr_square, which is the zz correlation square for all
 * configurations symmetric about spin chain center.
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
 * Initialize zz_all_corr_square in model_data_. The outer index is for J; middle index for for all possible
 * configurations, and the inner index is for realization
 */
void DisorderModelTransition::ZZ_all_corr_square_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N;
    const int num_realization = parameters.generic.num_realizations;
    const int size = parameters.generic.size;

    model_data_.zz_all_corr_square.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.zz_all_corr_square[i].resize(size/2);

        for (int j=0; j<size/2;j++){
            model_data_.zz_all_corr_square[i][j].resize(num_realization);

            for (int k=0; k<num_realization; k++){
                model_data_.zz_all_corr_square[i][j][k] = 0;
            }
        }
    }

    op_diag_ = true;
    op_dial_evec_keep_ = true;
}

/*
 * Compute zz_all_corr_square given realization index
 */
void DisorderModelTransition::ZZ_all_corr_square_compute_(AllPara const & parameters,
                                                          const EvolOP* model,
                                                          const DisorderLocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    const int dim = model -> Get_Dim();
    const int size = parameters.generic.size;

    if (local_info.evec_type_real){
        // Vector for eigenvectors in basic binary basis
        vector<vector<double> > evec_basic(dim);
        for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, local_info.evec_real, evec_basic);

        for (int i=0; i<size/2;i++){
            // For a given configuration
            const int right_up = 1 << i; // right reference spin
            const int left_up = 1 << (size - i - 1); // left reference spin
            // Record values for all eigenstate
            double zz_square = 0;

            for (int j=0; j<evec_basic.size();j++){
                // Record values for each eigenstate
                double left_temp = 0;
                double right_temp = 0;
                double left_right_temp = 0;

                for (int k=0; k<evec_basic[j].size(); k++){
                    int right_spin = -1;
                    int left_spin = -1;
                    if ( (k & right_up) == right_up) right_spin = 1;
                    if ( (k & left_up) == left_up ) left_spin = 1;

                    left_temp += left_spin * evec_basic[j][k] * evec_basic[j][k];
                    right_temp += right_spin * evec_basic[j][k] * evec_basic[j][k];
                    left_right_temp += left_spin * right_spin * evec_basic[j][k] * evec_basic[j][k];
                }

                zz_square += (left_right_temp - left_temp*right_temp ) * (left_right_temp - left_temp*right_temp );
                if (debug){
                    cout << "Realization " << local_info.realization_index << "  Configuration: " << i << endl;
                    cout << " eigenstate " << j << ":" << endl;
                    cout << "Average left spin: " << left_temp << endl;
                    cout << "Average right spin: " << right_temp << endl;
                    cout << "Average left X right: " << left_right_temp << endl;
                }
            }

            model_data_.zz_all_corr_square[local_info.J_index][i][local_info.realization_index]
                    += zz_square / double(dim);

            if (debug){
                cout << "Realization " << local_info.realization_index << "  Configuration: " << i << endl;
                cout << "Average correlation square: " << zz_square / double(dim) << endl;
                cout << endl;
            }
        }
    }
    else{
        // Vector for eigenvectors in basic binary basis
        vector<vector<complex<double> > > evec_basic(dim);
        for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, local_info.evec_complex, evec_basic);

        for (int i=0; i<size/2;i++){
            // For a given configuration
            const int right_up = 1 << i; // right reference spin
            const int left_up = 1 << (size - i - 1); // left reference spin
            // Record values for all eigenstate
            double zz_square = 0;

            for (int j=0; j<evec_basic.size();j++){
                // Record values for each eigenstate
                double left_temp = 0;
                double right_temp = 0;
                double left_right_temp = 0;

                for (int k=0; k<evec_basic[j].size(); k++){
                    int right_spin = -1;
                    int left_spin = -1;
                    if ( (k & right_up) == right_up) right_spin = 1;
                    if ( (k & left_up) == left_up ) left_spin = 1;

                    left_temp += left_spin * norm(evec_basic[j][k]);
                    right_temp += right_spin * norm(evec_basic[j][k]);
                    left_right_temp += left_spin * right_spin * norm(evec_basic[j][k]);
                }

                zz_square += (left_right_temp - left_temp*right_temp ) * (left_right_temp - left_temp*right_temp );
                if (debug){
                    cout << "Realization " << local_info.realization_index << "  Configuration: " << i << endl;
                    cout << " eigenstate " << j << ":" << endl;
                    cout << "Average left spin: " << left_temp << endl;
                    cout << "Average right spin: " << right_temp << endl;
                    cout << "Average left X right: " << left_right_temp << endl;
                }
            }

            model_data_.zz_all_corr_square[local_info.J_index][i][local_info.realization_index]
                    += zz_square / double(dim);

            if (debug){
                cout << "Realization " << local_info.realization_index << "  Configuration: " << i << endl;
                cout << "Average correlation square: " << zz_square / double(dim) << endl;
                cout << endl;
            }
        }
    }
}

/*
 * Output averages of mean end-to-end zz correlation squares (over all eigenstates for each realization) and their
 * standard deviations among all realizations for each configuration. For each J, a sperate file is outputted.
 */
void DisorderModelTransition::ZZ_all_corr_square_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.zz_all_corr_square.size() != J_N){
        cout << "Different number of J for zz all correlation square." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.zz_all_corr_square.size() << endl;
        abort();
    }

    for (int i=0; i<J_N; i++){
        if (model_data_.zz_all_corr_square[i].size() != size/2){
            cout << "Different number of configurations at " << i <<"th J for zz all correlation square." << endl;
            cout << "Expected Number: " << size/2 << endl;
            cout << "Actual Number: " << model_data_.zz_all_corr_square[i].size() << endl;
            abort();
        }

        for (int j=0; j<size/2;j++){
            if (model_data_.zz_all_corr_square[i][j].size() != num_realizations){
                cout << "Different number of realizations at " << i <<"th J and " << j << "th configuration"
                << " for zz all correlation square." << endl;
                cout << "Expected Number: " << num_realizations << endl;
                cout << "Actual Number: " << model_data_.zz_all_corr_square[i][j].size() << endl;
                abort();
            }
        }
    }

    stringstream prefix, suffix;
    prefix << name << ",size=" << size << ",Run=" << num_realizations << ",J=";
    suffix << ",zz_all_corr_square";
    if (version > 0) suffix <<",v" << version;
    suffix << ".txt";

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        stringstream filename;

        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        filename << prefix.str() << J << suffix.str();

        if (output) cout << filename.str() <<endl;

        ofstream fout( filename.str().c_str() );

        for (int j=0; j<model_data_.zz_all_corr_square[i].size(); j++){

            generic_mean_sd(model_data_.zz_all_corr_square[i][j], mean, sd);
            fout << setw(10) << j << setw(width) << mean << setw(width) << sd << endl;
        }
    }
}


