//
// Created by Liangsheng Zhang on 2/16/16.
//

/*
 * This file contains functions related to zz_full_time_corr, which is the full expression of end-to-end zz
 * time four point correlation used in studying floquet model transition from thermalization to localization.
 * The initial system is given a density matrix of identity matrix.
 *
 * Compared to the zz_time_corr, it has two extra terms:
 *     1/Z^2 sum_{n,m} ( |<n|sigma^z_0|m>|^2 |<n|sigma^z_L|m>|^2 + <n|sigma^z_0|m>^2 <m|sigma^z_L|n>^2 )
 * Note this term is real.
 */

#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include "disorder_model_transition.h"
#include "methods.h"
#include "generic_func.h"

using namespace std;

// Compute the extra two terms. Directly output the final result of the two terms
double zz_time_full_corr_extra_terms(const vector< vector<double> >&, bool );
double zz_time_full_corr_extra_terms(const vector< vector< complex<double> > >&, bool );


/*
 * Initialize zz_time_full_corr in model_data_. The outer index is for J, and the inner index is for realization
 */
void DisorderModelTransition::ZZ_time_full_corr_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.zz_time_full_corr.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.zz_time_full_corr[i].resize(num_realization);

        for (int j=0; j<num_realization;j++) model_data_.zz_time_full_corr[i][j] = 0;
    }

    op_diag_ = true;
    op_dial_evec_keep_ = true;
}

/*
 * Compute zz_full_time_corr given J index and realization index
 */
void DisorderModelTransition::ZZ_time_full_corr_compute_(AllPara const & parameters,
                                                    const EvolOP* model,
                                                    const DisorderLocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    const int dim = model -> Get_Dim();
    const int half_dim = dim / 2; // Below which leftmost spin is -1

    // Record values for all eigenstate
    double zz = 0;
    double z_left = 0;
    double z_right = 0;
    double z_left_ave_right_ave = 0;
    double extra_terms = 0;

    if (local_info.evec_type_real){
        // Vector for eigenvectors in basic binary basis
        vector<vector<double> > evec_basic(dim);
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

            zz += left_right_temp * left_right_temp;
            z_left += left_temp * left_temp;
            z_right += right_temp * right_temp;
            z_left_ave_right_ave += left_temp * right_temp;

            if (debug){
                cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
                cout << "Average left spin: " << left_temp << endl;
                cout << "Average right spin: " << right_temp << endl;
                cout << "Average left X right: " << left_right_temp << endl;
            }
        }

        extra_terms = zz_time_full_corr_extra_terms(evec_basic, debug);
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

            zz += left_right_temp * left_right_temp;
            z_left += left_temp * left_temp;
            z_right += right_temp * right_temp;
            z_left_ave_right_ave += left_temp * right_temp;

            if (debug){
                cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
                cout << "Average left spin: " << left_temp << endl;
                cout << "Average right spin: " << right_temp << endl;
                cout << "Average left X right: " << left_right_temp << endl;
            }
        }

        extra_terms = zz_time_full_corr_extra_terms(evec_basic, debug);

    }


    zz /= double(dim);
    z_left /= double(dim);
    z_right /= double(dim);
    z_left_ave_right_ave /= double(dim);

    double zz_time_full_corr = zz - z_left * z_right - z_left_ave_right_ave * z_left_ave_right_ave
                               - extra_terms;

    model_data_.zz_time_full_corr[local_info.J_index][local_info.realization_index] += zz_time_full_corr;
    if (debug){
        cout << "Realization " << local_info.realization_index << " average time correlation: "
        << zz_time_full_corr << endl;
        cout << endl;
    }
}

/*
 * Output averages of mean end-to-end zz time four-point correlation squares and their standard deviations among all
 * realizations for each J
 */
void DisorderModelTransition::ZZ_time_full_corr_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.zz_time_full_corr.size() != J_N){
        cout << "Not enough number of J for zz time full correlation." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.zz_time_full_corr.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.zz_time_full_corr[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for zz time full correlation."
                 << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.zz_time_full_corr[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",zz_time_full_corr";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        generic_mean_sd(model_data_.zz_time_full_corr[i], mean, sd);
        fout << setw(10) << J << setw(width) << mean << setw(width) << sd << endl;
    }
}

double zz_time_full_corr_extra_terms(const vector< vector<double> >& evec_basic, bool debug){
    const int dim = evec_basic.size();
    const int half_dim = dim/2;
    double extra_terms = 0;

    for (int i=0; i<evec_basic.size(); i++){
        for (int j=0; j<evec_basic.size(); j++){
            // Record values for each pair of eigenstate
            double left_temp = 0;
            double right_temp = 0;

            if(evec_basic[i].size() != dim){
                cout << "Incompatible dimension for eigenstate " << i << " in extra terms for zz time full corr"
                     << endl;
                cout << "Obtained: " << evec_basic[i].size() << endl;
                cout << "Expected: " << dim << endl;
                abort();
            }

            for (int k=0; k<evec_basic[i].size(); k++){
                int right_spin = -1;
                int left_spin = 1;
                if (k & 1 == 1) right_spin = 1;
                if (k < half_dim ) left_spin = -1;

                left_temp += left_spin * evec_basic[i][k] * evec_basic[j][k];
                right_temp += right_spin * evec_basic[i][k] * evec_basic[j][k];
            }

            // The two extra terms are equal in the real case
            extra_terms += 2 * left_temp * left_temp * right_temp * right_temp;

            if (debug){
                cout << "eigenstate " << i << " and " << j << ":" << endl;
                cout << "Left correlation: " << left_temp << endl;
                cout << "Right correlation: " << right_temp << endl;
            }
        }
    }

    return extra_terms / double(dim*dim);
}

double zz_time_full_corr_extra_terms(const vector< vector< complex<double> > >& evec_basic, bool debug){
    const int dim = evec_basic.size();
    const int half_dim = dim/2;
    double extra_terms = 0;

    for (int i=0; i<evec_basic.size(); i++){
        for (int j=0; j<evec_basic.size(); j++){
            // Record values for each pair of eigenstate
            complex<double> left_temp = 0;
            complex<double> right_temp = 0;

            if(evec_basic[i].size() != dim){
                cout << "Incompatible dimension for eigenstate " << i << " in extra terms for zz time full corr"
                << endl;
                cout << "Obtained: " << evec_basic[i].size() << endl;
                cout << "Expected: " << dim << endl;
                abort();
            }

            for (int k=0; k<evec_basic[i].size(); k++){
                int right_spin = -1;
                int left_spin = 1;
                if (k & 1 == 1) right_spin = 1;
                if (k < half_dim ) left_spin = -1;

                left_temp += left_spin * conj(evec_basic[i][k]) * evec_basic[j][k];
                right_temp += right_spin * evec_basic[i][k] * conj(evec_basic[j][k]);
            }

            complex<double> sqr = left_temp*left_temp*right_temp*right_temp;
            extra_terms += norm(left_temp * right_temp) + sqr.real();

            if (debug){
                cout << "eigenstate " << i << " and " << j << ":" << endl;
                cout << "Left correlation: " << left_temp << endl;
                cout << "Right correlation: " << right_temp << endl;
            }
        }
    }

    return extra_terms / double(dim*dim);
}





