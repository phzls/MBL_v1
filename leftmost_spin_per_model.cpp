//
// Created by Liangsheng Zhang on 6/12/15.
//

#include <cmath>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "evol_data.h"
#include "methods.h"
#include "screen_output.h"

using namespace std;
using namespace Eigen;

/*
 * This file contains functions related to leftmost_spin_per_model. It is the average value of the leftmost spin
 * at the spin chain. When outputted, the infinite time value computed from diagonal ensemble is subtracted.
 */

/*
 * Initialize the leftmost spin. The outer index is model, middle index is time and the inner index is realization.
 * Initialize the infinite time leftmost spin. The outter index is model and the inner index is realization.
 */
void EvolData::Leftmost_Spin_Per_Model_Init_(const AllPara& parameters) {
    const int num_realizations = 1; // Since no random number in initial state, for now only 1 realization
    const int time_step = parameters.evolution.time_step;

    leftmost_spin_per_model_.resize(parameters.evolution.model_num);
    leftmost_spin_infinite_time_per_model_.resize(parameters.evolution.model_num);

    for (int i=0; i<parameters.evolution.model_num; i++){
        leftmost_spin_per_model_[i].resize(time_step);
        leftmost_spin_infinite_time_per_model_[i].resize(num_realizations);

        for (int j=0; j<num_realizations; j++) leftmost_spin_infinite_time_per_model_[i][j] = 0;

        for (int j=0; j<time_step; j++){
            leftmost_spin_per_model_[i][j].resize(num_realizations);

            for (int k=0; k<num_realizations; k++){
                leftmost_spin_per_model_[i][j][k] = 0;
            }
        }
    }
}

/*
 * Use information from initial condition to find infinite time leftmost spin
 */
void EvolData::Leftmost_Spin_Per_Model_Evol_Init_(const InitEvolData& init_evol_data,
                                                  const InitEvolInfo& init_evol_info) {
    const int model = init_evol_info.model;
    const int realization = init_evol_info.realization;

    leftmost_spin_infinite_time_per_model_[model][realization] = init_evol_data.infinite_time_leftmost_spin;
}

/*
 * Compute the leftmost spin at a given time step and realization, given a state density matrix
 */
void EvolData::Leftmost_Spin_Per_Model_Cal_C_(const MatrixXcd& density_matrix, const StepInfo& info){
    const int realization = info.realization;
    const int time = info.time;
    const int model = info.model;
    const string basis_type = info.basis_type;
    const int total_rank = density_matrix.rows();
    const int rest_rank = total_rank >> 1;

    if (realization == 1){
        if ((basis_type == "Binary") || (basis_type == "binary")){
            // Check whether density_matrix is Hermitian
            if (density_matrix.rows() != density_matrix.cols()){
                cout << "Density matrix passed in entropy_per_model_cal_C is not square." << endl;
                cout << "Rows: " << density_matrix.rows() << endl;
                cout << "Cols: " << density_matrix.cols() << endl;
                abort();
            }

            for (int i=0; i< density_matrix.rows(); i++){
                for (int j=i; j<density_matrix.rows();j++){
                    if ( norm( density_matrix(i,j) - conj(density_matrix(j,i)) ) > info.delta ){
                        cout << "Density matrix is not Hermitian at (" << i << "," << j << ")." << endl;
                        cout << "At (" << i << "," << j << "): " << density_matrix(i,j) << endl;
                        cout << "At (" << j << "," << i << "): " <<  density_matrix(j,i) << endl;
                        abort();
                    }
                }
            }

            double leftmost_spin = 0;

            for (int i=0; i<rest_rank;i++){
                leftmost_spin -= real(density_matrix(i,i));
            }
            for (int i=rest_rank; i<total_rank; i++){
                leftmost_spin += real(density_matrix(i,i));
            }

            leftmost_spin_per_model_[model][time][realization] = leftmost_spin;

            if (info.debug){
                cout << "Leftmost spin per model:" << endl;
                cout << leftmost_spin_per_model_[model][time][realization] << endl;
                cout << endl;
            }
        }
        else if ((basis_type == "Evec") || (basis_type == "evec")){
            cout << "Density matrix in evolution eigenvector basis is not "
            << "used for leftmost spin computation." << endl;
        }
        else{
            cout << "Unknown basis type " << basis_type << endl;
            abort();
        }
    }
}

/*
 * Compute the leftmost spin at a given time step and realization, given a state vector
 */
void EvolData::Leftmost_Spin_Per_Model_Cal_(const VectorXcd& state, const StepInfo& info){
    cout << "Not implemented yet." << endl;
    abort();
}

/*
 * Output leftmost spin minus corresponding infinite time value per model. The spin difference will
 * be averaged over different realizations, and the sample variance will be computed, where the sample
 * is taken at a fixed time step with different realizations. The time step, average and the standard
 * deviation computed from the sample variance will be computed. Note this standard deviation is not
 * the standard deviation of average yet. The name taken in will become part of the output file name.
 */

void EvolData::Leftmost_Spin_Per_Model_Out_(const AllPara& parameters, const string& name){
    const int time_step = parameters.evolution.time_step;
    const int model_num = parameters.evolution.model_num;
    const int num_realizations = parameters.generic.num_realizations;
    const int jump = parameters.evolution.jump;
    const bool output = parameters.output.filename_output;
    const bool log_time = parameters.evolution.log_time;
    const bool markov_jump = parameters.evolution.markov_jump;
    const int markov_time_jump = parameters.evolution.markov_time_jump;
    const int version = parameters.generic.version;

    if (leftmost_spin_per_model_.size() != model_num){
        cout << "Not enough models are computed for leftmost spin." << endl;
        cout << "Expected: " << model_num << endl;
        cout << "Computed: " << leftmost_spin_per_model_.size() << endl;
    }

    if (leftmost_spin_infinite_time_per_model_.size() != model_num){
        cout << "Not enough models are computed for infinite time leftmost spin." << endl;
        cout << "Expected: " << model_num << endl;
        cout << "Computed: " << leftmost_spin_per_model_.size() << endl;
    }

    for (int i=0; i<model_num;i++){
        if (leftmost_spin_per_model_[i].size() != time_step){
            cout << "Not enough time steps are computed for leftmost spin at " << i << "th model." << endl;
            cout << "Expected: " << time_step << endl;
            cout << "Computed: " << leftmost_spin_per_model_[i].size() << endl;
        }

        if (leftmost_spin_infinite_time_per_model_[i].size() != num_realizations){
            cout << "Not enough realizations are computed for infinite time leftmost spin at "
                 << i << "th model." << endl;
            cout << "Expected: " << num_realizations << endl;
            cout << "Computed: " << leftmost_spin_infinite_time_per_model_[i].size() << endl;
        }
    }

    for (int i=0; i<model_num;i++){
        for(int j=0; j<time_step; j++){
            if (num_realizations != leftmost_spin_per_model_[i][j].size()){
                cout << "Entropy is not computed with enough realizations at " << i
                << "th model for time " << i << endl;
                cout << "Expect: " << num_realizations << endl;
                cout << "Computed: " << leftmost_spin_per_model_[i][j].size() << endl;
                abort();
            }
        }
    }

    stringstream filename;
    filename << name <<",Run=" << num_realizations << ",Total=" << time_step;

    if (log_time) filename <<",log_time";

    if (markov_jump) filename <<",markov_jump=" << markov_time_jump;

    filename << ",jump=" << jump << ",leftmost_spin_difference_per_model";

    if (parameters.evolution.sample_detail) filename << "_sample_detail";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    // Subtract the infinite time value
    for (int i=0; i<model_num;i++){
        for (int j=0; j<num_realizations; j++){
            double inf_spin = leftmost_spin_infinite_time_per_model_[i][j];

            for (int k=0; k<time_step; k++){
                leftmost_spin_per_model_[i][j][k] -= inf_spin;
            }
        }
    }

    General_Output_(parameters, leftmost_spin_per_model_, filename.str());
}


