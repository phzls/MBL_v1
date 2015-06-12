//
// Created by Liangsheng Zhang on 6/1/15.
//

#include <cmath>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "evol_data.h"
#include "methods.h"
#include "generic_func.h"
#include "screen_output.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to entropy_per_model. It is the entropy of the
 ** left part of a spin chain model whose local dimension is 2. The entropy is computed in bits.
 **/

/*
 * Initialize the entropy. The outer index is model, middle index is time and the inner index is realization.
 */
void EvolData::Entropy_Per_Model_Init_(const AllPara& parameters){
    const int num_realizations = parameters.generic.num_realizations;
    const int time_step = parameters.evolution.time_step;

    entropy_per_model_.resize(parameters.evolution.model_num);

    for (int i=0; i<parameters.evolution.model_num; i++){
        entropy_per_model_[i].resize(time_step);

        for (int j=0; j<time_step; j++){
            entropy_per_model_[i][j].resize(num_realizations);

            for (int k=0; k<num_realizations; k++){
                entropy_per_model_[i][j][k] = 0;
            }
        }
    }
}

/*
 * Use information from initial states for evolution. So far it is not used.
 */
void EvolData::Entropy_Per_Model_Evol_Init_(const InitEvolData& init_evol_data,
                                            const InitEvolInfo& init_evol_info) {
    if (init_evol_info.debug){
        cout << "Init_evol_data is not used in entropy_per_model." << endl;
    }
}

/*
 * Compute the entropy at a given time step and realization, given a state vector
 */
void EvolData::Entropy_Per_Model_Cal_(const VectorXcd& state_basic, const StepInfo& info){
    const int realization = info.realization;
    const int time = info.time;
    const int model = info.model;
    const int left_size = info.left_size;
    const string basis_type = info.basis_type;

    if ((basis_type == "Binary") || (basis_type == "binary")){
        MatrixXcd reduced_density; // Reduced density matrix for the left part
        reduced_density_left_2(state_basic, size_, left_size, reduced_density);

        if (info.debug){
            cout << "Reduced density matrix:" << endl;
            complex_matrix_write(reduced_density);
            cout << endl;
        }

        SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix

        density_eigen.compute(reduced_density, false); // Eigenvectors not computed

        entropy_per_model_[model][time][realization] = 0;

        for (int i=0; i<density_eigen.eigenvalues().rows();i++){
            double eval = density_eigen.eigenvalues()(i);

            if (abs(eval)>1.0e-15)
            {
                if (eval<0){
                    cout << "Density matrix has significant negative eigenvalues." << endl;
                    cout << eval << endl;
                    abort();
                }
                entropy_per_model_[model][time][realization] += -eval*log2(eval);
            }
        }

        if (info.debug){
            cout << "Entropy per model:" << endl;
            cout << entropy_per_model_[model][time][realization] << endl;
            cout << endl;
        }
    }
    else if ((basis_type == "Evec") || (basis_type == "evec")){
        if (info.debug){
            cout << "State vector in evolution eigenvector basis"
                 << " is not used for half-chain entropy computation." << endl;
        }
    }
    else{
        cout << "Unknown basis type " << basis_type << endl;
        abort();
    }
}


/*
 * Compute the entropy at a given time step and realization, given a complex density matrix.
 * Entropy is for the whole system.
 */
void EvolData::Entropy_Per_Model_Cal_C_(const MatrixXcd& density_matrix, const StepInfo& info){
    const int realization = info.realization;
    const int time = info.time;
    const int model = info.model;
    const string basis_type = info.basis_type;

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

        SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix

        density_eigen.compute(density_matrix, false); // Eigenvectors not computed

        entropy_per_model_[model][time][realization] = 0;

        for (int i=0; i<density_eigen.eigenvalues().rows();i++){
            double eval = density_eigen.eigenvalues()(i);

            if (abs(eval)>1.0e-8)
            {
                if ( eval*log2(eval) != eval*log2(eval) ){
                    cout << "Significant negative eigenvalues of density matrix." << endl;
                    cout << "eval: " << eval << endl;
                    cout << "Time: " << time << "  Realization: " << realization << endl;
                    abort();
                }
                entropy_per_model_[model][time][realization] += -eval*log2(eval);
            }
        }

        if (info.debug){
            cout << "Entropy per model:" << endl;
            cout << entropy_per_model_[model][time][realization] << endl;
            cout << endl;
        }
    }
    else if ((basis_type == "Evec") || (basis_type == "evec")){
        cout << "Density matrix in evolution eigenvector basis is not "
             << "used for full-chain entropy computation." << endl;
    }
    else{
        cout << "Unknown basis type " << basis_type << endl;
        abort();
    }
}

/*
 * Output entropy per model. The entroy will be averaged over different realizations, and the
 * sample variance will be computed, where the sample is taken at a fixed time step with different
 * realizations. The time step, average and the standard deviation computed from the sample
 * variance will be computed. Note this standard deviation is not the standard deviation of
 * average yet. The name taken in will become part of the output file name.
 */

void EvolData::Entropy_Per_Model_Out_(const AllPara& parameters, const string& name){
    const int time_step = parameters.evolution.time_step;
    const int model_num = parameters.evolution.model_num;
    const int num_realizations = parameters.generic.num_realizations;
    const int jump = parameters.evolution.jump;
    const bool output = parameters.output.filename_output;
    const int width = parameters.output.width;
    const double step_size = parameters.evolution.step_size;
    const bool log_time = parameters.evolution.log_time;
    const int log_time_jump = parameters.evolution.log_time_jump;
    const bool markov_jump = parameters.evolution.markov_jump;
    const int markov_time_jump = parameters.evolution.markov_time_jump;
    const int version = parameters.generic.version;

    if (entropy_per_model_.size() != model_num){
        cout << "Not enough models are computed for entropy." << endl;
        cout << "Expected: " << model_num << endl;
        cout << "Computed: " << entropy_per_model_.size() << endl;
    }

    for (int i=0; i<model_num;i++){
        if (entropy_per_model_[i].size() != time_step){
            cout << "Not enough time steps are computed for entropy at " << i << "th model." << endl;
            cout << "Expected: " << time_step << endl;
            cout << "Computed: " << entropy_per_model_[i].size() << endl;
        }
    }

    for (int i=0; i<model_num;i++){
        for(int j=0; j<time_step; j++){
            if (num_realizations != entropy_per_model_[i][j].size()){
                cout << "Entropy is not computed with enough realizations at " << i
                     << "th model for time " << i << endl;
                cout << "Expect: " << num_realizations << endl;
                cout << "Computed: " << entropy_per_model_[i][j].size() << endl;
                abort();
            }
        }
    }

    stringstream filename;
    filename << name <<",Run=" << num_realizations << ",Total=" << time_step;

    if (log_time) filename <<",log_time";

    if (markov_jump) filename <<",markov_jump=" << markov_time_jump;

    filename << ",jump=" << jump << ",entropy_per_model";

    if (parameters.evolution.sample_detail) filename << "_sample_detail";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int t=0; t<time_step; t++){
        double time = t * step_size * jump; // An overflow still happens for entropy

        if (markov_jump) time *= markov_time_jump;

        if (log_time){
            long long int power = pow(log_time_jump,t);
            time = step_size * power;
        }

        if (parameters.evolution.sample_detail && (entropy_per_model_.size() == 1)){
            fout << setw(10) << time;

            for (int i=0; i<entropy_per_model_[0][t].size(); i++) fout << setw(width) << entropy_per_model_[0][t][i];
            fout << endl;
        }
        else if (parameters.evolution.sample_detail){
            fout << setw(10) << time;
            for (int i=0; i < entropy_per_model_.size(); i++){
                double mean, sd;
                generic_mean_sd(entropy_per_model_[i][t], mean, sd);
                // For now errors are not outputted
                fout << setw(width) << mean;
            }
            fout << endl;
        }
        else{
            const int model_num = entropy_per_model_.size();
            vector<double> mean(model_num);
            vector<double> sd(model_num);

            for (int i=0; i<model_num; i++){
                generic_mean_sd(entropy_per_model_[i][t], mean[i], sd[i]);
            }

            double final_mean, final_sd;
            generic_mean_sd(mean, final_mean, final_sd);

            fout << setw(10) << time << setw(width) << final_mean << setw(width) << final_sd << endl;
        }
    }
}

