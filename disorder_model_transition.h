//
// Created by Liangsheng Zhang on 4/15/15.
//

#ifndef MBL_V1_DISORDER_MODEL_TRANSITION_H
#define MBL_V1_DISORDER_MODEL_TRANSITION_H

#include <iostream>
#include <vector>
#include <map>
#include "parameters.h"
#include "evol_op.h"

using namespace std;

/*
 * The data class used in flo_transition. It contains all possible outputs.
 */

struct DisorderModelData
{
    // End-to-end sigma_z-sigma_z correlation square averaged over all eigenstates and realizations
    vector< vector<double> > zz_corr_square;

    // Entropy varaince for all eigenstates
    vector< vector<double> > ent_var;

    // Entropy variance for the eigenstate with smallest phase magnitude among different realizations
    vector< vector<double> > ent_smallest_var;

    // End-to-end sigma_z-sigma_z time four-point correlation
    vector< vector<double> > zz_time_corr;

    // Various components of zz time correlation
    vector< vector<vector<double> > > zz_time_corr_component;

    // zz correlation at different distances with configuration symmetric w.r.t the center
    vector< vector<vector<double> > > zz_all_corr_square;

    // End-to-end sigma_z-sigma_z correlation square with its logarithm averaged over all eigenstates and realizations
    vector< vector<double> > log_zz_corr_square;

    // log zz correlation square at different distances with configuration symmetric w.r.t the center
    vector< vector<vector<double> > > log_zz_all_corr_square;

    // zz time correlation at different distances with configuration symmetric w.r.t the center
    vector< vector<vector<double> > > zz_all_time_corr;

    // Mean value of entropy for each realization
    vector< vector<double> > ent_mean;

    // Number of eigenstates, assuming it's the same across the calculation
    int model_dim;

    // Level statistics average
    vector< vector<double> > level_stat_ave;
};

/*
 * Local information that helps computation
 */
struct DisorderLocalInfo
{
    int J_index; // The index of J
    int realization_index; // The index of realization

    bool evec_type_real; // Whether real eigenvectors are recorded
    bool eval_type_real; // Whether real eigenvalues are recorded

    vector<MatrixXcd> evec_complex; // Eigenvectors
    vector<MatrixXd> evec_real; // Eigenvectors

    vector<VectorXcd> eval_complex; // Eigenvalues
    vector<VectorXd> eval_real; // Eigenvalues
};

class DisorderModelTransition;

// Pointer to functions initializing data
typedef  void (DisorderModelTransition::*Flo_init)(const AllPara&);

// Pointer to functions studying properties of a given floquet system
typedef void (DisorderModelTransition::*Flo_func)(const AllPara&, const EvolOP*, const DisorderLocalInfo&);

// Pointer to functions outputting data. The string gives the first part of filename
typedef void (DisorderModelTransition::*Flo_out)(const AllPara&, const string&);

/*
 * The class that contains data and functions which handle it
 */
class DisorderModelTransition
{
private:
    // Map that determines which functions are computed
    map<string, bool> flo_func_bool_map_;

    // Map that links strings to initialization functions
    map<string, Flo_init> flo_init_map_;

    // Map that links strings to computation functions
    map<string, Flo_func> flo_func_map_;

    // Map that links strings to output functions
    map<string, Flo_out> flo_out_map_;

    // Whether eigenvectors can be real
    bool evec_type_real_;

    // Whether eigenvalues can be real
    bool eval_type_real_;

    // Whether keep evolution operator
    bool op_keep_;

    // Whether diagonalize evolution operator
    bool op_diag_;

    // Whether keep eigenvectors when diagonalization
    bool op_dial_evec_keep_;

    // Initialize all maps
    void map_initialize_(const AllPara&);

    DisorderModelData model_data_;

    // For end to end zz correlation square
    void ZZ_corr_square_init_(const AllPara&);
    void ZZ_corr_square_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void ZZ_corr_square_out_(const AllPara&, const string&);

    // Output end to end zz correlation square for all samples.
    // It calls same functions as zz correlation square for iniitalization and computation
    // if they are not called
    void ZZ_corr_square_all_sample_init_(const AllPara&);
    void ZZ_corr_square_all_sample_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void ZZ_corr_square_all_sample_out_(const AllPara&, const string&);

    // For variance of entropy among eigenstates
    void Ent_var_init_(const AllPara&);
    void Ent_var_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void Ent_var_out_(const AllPara&, const string&);

    // For end to end zz time four-point correlation
    void ZZ_time_corr_init_(const AllPara&);
    void ZZ_time_corr_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void ZZ_time_corr_out_(const AllPara&, const string&);

    // Output end to end zz four-point correlation for all samples.
    // It calls same functions as zz four-point correlation for iniitalization and computation
    // if they are not called
    void ZZ_time_corr_all_sample_init_(const AllPara&);
    void ZZ_time_corr_all_sample_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void ZZ_time_corr_all_sample_out_(const AllPara&, const string&);

    // For end to end zz time four-point correlation components
    void ZZ_time_corr_component_init_(const AllPara&);
    void ZZ_time_corr_component_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void ZZ_time_corr_component_out_(const AllPara&, const string&);

    // For zz correlation square at all distances
    void ZZ_all_corr_square_init_(const AllPara&);
    void ZZ_all_corr_square_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void ZZ_all_corr_square_out_(const AllPara&, const string&);

    // For log end to end zz correlation square, used in MBL phase
    void Log_ZZ_corr_square_init_(const AllPara&);
    void Log_ZZ_corr_square_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void Log_ZZ_corr_square_out_(const AllPara&, const string&);

    // For log zz correlation square at all distances
    void Log_ZZ_all_corr_square_init_(const AllPara&);
    void Log_ZZ_all_corr_square_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void Log_ZZ_all_corr_square_out_(const AllPara&, const string&);

    // For zz time correlation at all distances
    void ZZ_all_time_corr_init_(const AllPara&);
    void ZZ_all_time_corr_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void ZZ_all_time_corr_out_(const AllPara&, const string&);

    // For "variance" of entropy among eigenstates, where the mean used in variance calculation is
    // the mean from all samples of all realizations.
    void Ent_var_all_mean_init_(const AllPara&);
    void Ent_var_all_mean_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void Ent_var_all_mean_out_(const AllPara&, const string&);
    bool ent_mean_keep_; // Whether keep mean values of entropy for each realization

    // For average (consecutive) level statistics for each realization from floquet system
    void Flo_level_stats_ave_init_(const AllPara&);
    void Flo_level_stats_ave_compute_(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void Flo_level_stats_ave_out_(const AllPara&, const string&);

public:
    DisorderModelTransition(const AllPara& parameters) : evec_type_real_(true), eval_type_real_(true),
                                                         op_keep_(false), op_diag_(false), op_dial_evec_keep_(false),
                                                         ent_mean_keep_(false)
    {
        map_initialize_(parameters);

        if (op_dial_evec_keep_) op_diag_ = true;
    }

    void Compute(const AllPara&, const EvolOP*, const DisorderLocalInfo&);
    void Output(const AllPara&, const string&);

    bool Evec_Real() const {return evec_type_real_;} // Whether evec can be real
    bool Eval_Real() const {return eval_type_real_;} // Whether eval can be real

    bool Op_Keep() const {return op_keep_;} // Whether evolution operator will be kept
    bool Op_Diag() const {return op_diag_;} // Whether evolution operator will be diagonalized
    bool Op_Evec_Keep() const {return op_dial_evec_keep_;} // Whether evec will be computed

    ~DisorderModelTransition() {};
};

#endif //MBL_V1_DISORDER_MODEL_TRANSITION_H
