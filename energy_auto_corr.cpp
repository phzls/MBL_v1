//
// Created by Liangsheng Zhang on 9/25/15.
//

/*
 * This file implements functions related to computing energy autocorrelation
 */

#include <cmath>
#include <fstream>
#include <complex>
#include <iomanip>
#include <sstream>
#include "flo_op_auto_corr.h"
#include "energy_auto_corr_func.h"
#include "screen_output.h"
#include "generic_func.h"

using namespace std;

/*
 * Initialize the vector which holds the data. The innermost index is for time. The middle
 * index is for realizations. The outer index is for different choices of certain parameters
 */
void OpAutoCorr::Energy_corr_init_(const AllPara& parameters) {
    const int num_realizations = parameters.generic.num_realizations;
    const int time_pts = parameters.flo_op_auto_corr.time_pts;
    const int para_pts = parameters.flo_op_auto_corr.para_pts;

    op_corr_data_.ene_auto_corr.resize(para_pts);

    for(int i=0; i<para_pts;i++){
        op_corr_data_.ene_auto_corr[i].resize(num_realizations);
        for(int j=0; j<num_realizations;j++){
            op_corr_data_.ene_auto_corr[i][j].resize(time_pts);
            for(int k=0; k<num_realizations;k++){
                op_corr_data_.ene_auto_corr[i][j][k] = 0;
            }
        }
    }

    // For now not checking duplicates
    string type = "Ising_Random_Simp_Shift_Cos_Real_Tau_Floquet";
    op_corr_data_.ham_construct_map_[type] = Ham_for_Ising_Random_Simp_Shift_Cos_Real_Flo;
}

/*
 * Set up the Hamiltonian matrix
 */
void OpAutoCorr::Energy_corr_set_up_(const AllPara& parameters, const EvolOP* evol) {
    const bool debug = parameters.generic.debug;
    const string type = evol -> Type();
    (*op_corr_data_.ham_construct_map_[type])(evol, op_corr_data_.ham);

    if(debug){
        cout << "Constructed Hamiltonian: " << endl;
        complex_matrix_write(op_corr_data_.ham);
        cout << endl;
    }
}

/*
 * Compute autocorrelation at the given time point
 */
void OpAutoCorr::Energy_corr_compute_(const AllPara& parameters, const OpCorrLocalInfo& local_info) {
    const int realization_num = local_info.realization_num;
    const int time_step = local_info.time_step;
    const int para_num = local_info.para_num;
    const int dim = op_corr_data_.ham.rows();
    const bool debug = parameters.generic.debug;

    vector< complex<double> > eval_power(dim);
    int index = 0;
    for(int i=0; i< eval.size(); i++){
        for(int j=0; j<eval[i].rows(); j++){
            if(index >= dim){
                cout << "Too many eigenvalues." << endl;
                cout << "Expected number: " << dim << endl;
                abort();
            }

            eval_power[index] = pow( eval[i](j), double(time_step) );
            ++ index;
        }
    }

    if(index != dim){
        cout << "Incompatible eigenvalue number." << endl;
        cout << "Hilbert dimension: " << dim << endl;
        cout << "Number of eigenvalues obtained: " << index << endl;
        abort();
    }

    // Compute trace
    complex<double> trace = 0;
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            trace += norm(op_corr_data_.ham(i,j)) * eval_power[i] * conj(eval_power[j]);
        }
    }

    if(abs(imag(trace)) > 1.0e-8){
        cout << "Trace for energy autocorrelation should be real." << endl;
        cout << "Obtained trace: " << trace << endl;
        abort();
    }

    op_corr_data_.ene_auto_corr[para_num][realization_num][time_step] = real(trace);

    if(debug){
        cout << "At realization " << realization_num << " and time step " << time_step << endl;
        cout << "Eigenvalue power:" << endl;
        complex_write(eval_power);
        cout << endl << "Trace: " << real(trace) << endl;
    }
}

/*
 * Results output: the first column would be for different choices of parameters,
 * at each row, autocorrelations are listed according to time
 * Mean and SD are recorded in two different files.
 */
void OpAutoCorr::Energy_corr_out_(const AllPara& parameters, const OpCorrLocalPara& local_para,
                                  const string& name) {
    const int para_pts = parameters.flo_op_auto_corr.para_pts;
    const int time_pts = parameters.flo_op_auto_corr.time_pts;
    const int num_realizations = parameters.generic.num_realizations;
    const double min_para = local_para.Min_Para();
    const double max_para = local_para.Max_Para();
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (op_corr_data_.ene_auto_corr.size() != para_pts){
        cout << "Incompatible number of J for energy autocorrelation." << endl;
        cout << "Expected Number: " << para_pts << endl;
        cout << "Actual number: " << op_corr_data_.ene_auto_corr.size() << endl;
        abort();
    }

    for (int i=0; i< para_pts; i++){
        if (op_corr_data_.ene_auto_corr[i].size() != num_realizations){
            cout << "Incompatible number of realizations at " << i <<"th parameter for energy autocorrelation." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << op_corr_data_.ene_auto_corr[i].size() << endl;
            abort();
        }
    }

    for(int i=0; i<para_pts; i++){
        for(int j=0; j<num_realizations; j++){
            if(op_corr_data_.ene_auto_corr[i][j].size() != time_pts){
                cout << "Incompatible number of time points at " << i << "th parameter and " << j << "th realization"
                     << " for energy autocorrelation." << endl;
                cout << "Expected Number: " << time_pts << endl;
                cout << "Actual Number: " << op_corr_data_.ene_auto_corr[i][j].size() << endl;
            }
        }
    }

    stringstream base_filename;
    stringstream mean_filename; // For mean value
    stringstream sd_filename; // For SD
    base_filename << name << ",size=" << size << ",para_N=" << para_pts << ",para_min=" << min_para
             << ",para_max=" << max_para << ",Run=" << num_realizations << ",N_time=" << time_pts
             << "ene_auto_corr";
    mean_filename << base_filename << "_mean";
    sd_filename << base_filename << "_sd";

    if (version > 0) {
        mean_filename <<",v" << version;
        sd_filename << ",v" << version;
    }

    mean_filename << ".txt";
    sd_filename << ".txt";

    if (output){
        cout << "Mean value file:" << endl << mean_filename.str().c_str() << endl;
        cout << "SD value file:" << endl << sd_filename.str().c_str() << endl;
    }

    ofstream fout_mean( mean_filename.str().c_str() );
    ofstream fout_sd( sd_filename.str().c_str() );

    for (int i=0; i<para_pts;i++){
        double para;
        if (para_pts > 1) para = min_para + i * (max_para - min_para)/(para_pts-1);
        else para = min_para;

        fout_mean << setw(10) << para;
        fout_sd << setw(10) << para;

        // Store value at same parameter and time point for different realizations
        vector<double> val(num_realizations);

        for(int j=0; j<time_pts; j++){
            double mean,sd;

            for(int k=0; k<num_realizations;k++){
                val[k] = op_corr_data_.ene_auto_corr[i][k][j];
            }

            generic_mean_sd(val, mean, sd);

            fout_mean << setw(width) << mean;
            fout_sd << setw(width) << sd;
        }

        fout_mean << endl;
        fout_sd << endl;
    }
}


