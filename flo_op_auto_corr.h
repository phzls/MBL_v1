//
// Created by Liangsheng Zhang on 9/24/15.
//

#ifndef MBL_V1_FLO_OP_AUTO_CORR_H
#define MBL_V1_FLO_OP_AUTO_CORR_H

#include <iostream>
#include <vector>
#include <map>
#include "parameters.h"
#include "evol_op.h"

using namespace std;
using namespace Eigen;

// Pointer to the function that constructs the operator for autocorrelation
typedef void (*Op_construct)(const EvolOP*, MatrixXcd&);

/*
 * The class that contains operators for autocorrelation and corresponding
 * autocorrelations
 */
struct OpCorrData
{
    MatrixXcd ham; // Hamiltonian for energy autocorrelation

    // Vectors holding energy autocorrelation. Innermost index is time, middle for realizations
    // and outer index for number of changing parameters
    vector< vector< vector<double> > > ene_auto_corr;
    map<string, Op_construct> ham_construct_map_; // Map from name to Hamiltonian constructed

    OpCorrData(){ ham.resize(0,0);}

    ~OpCorrData() { ham.resize(0,0); }
};

/*
 * Local information that helps computation
 */
struct OpCorrLocalInfo
{
    int time_step; // Index for time
    int para_num; // Index for parameter loop
    int realization_num; // Index for realization
};

class OpAutoCorr;
class OpCorrLocalPara;

// Pointer to functions that initialize the data
typedef void (OpAutoCorr::*Op_init)(const AllPara&);

// Pointer to functions that set up the operator
typedef  void (OpAutoCorr::*Op_set_up)(const AllPara&, const EvolOP*);

// Pointer to functions studying properties of a given model
typedef void (OpAutoCorr::*Op_corr_cal)(const AllPara&, const OpCorrLocalInfo&);

// Pointer to functions outputting data. The string gives the first part of filename
typedef void (OpAutoCorr::*Op_corr_out)(const AllPara&, const OpCorrLocalPara&, const string&);

/*
 * The class that contains data and functions which handle it
 */
class OpAutoCorr
{
private:
    // Map that determines which functions are computed
    map<string, bool> op_auto_corr_map_;

    // Map that links strings to initialization functions
    map<string, Op_init> op_init_map_;

    // Map that links strings to setting up functions
    map<string, Op_set_up> op_set_up_map_;

    // Map that links strings to computation functions
    map<string, Op_corr_cal> op_corr_cal_map_;

    // Map that links strings to output functions
    map<string, Op_corr_out> op_corr_out_map_;

    OpCorrData op_corr_data_;

    // Compute energy autocorrelation
    void Energy_corr_init_(const AllPara&);
    void Energy_corr_set_up_(const AllPara&, const EvolOP*);
    void Energy_corr_compute_(const AllPara&, const OpCorrLocalInfo&);
    void Energy_corr_out_(const AllPara&, const OpCorrLocalPara&, const string&);

public:
    vector<VectorXcd> eval; // Used to record eigenvalues; for now complex is enough

    // Initialize all maps
    OpAutoCorr(const AllPara&);

    // Set up calculation
    void SetUp(const AllPara&, const EvolOP*);

    // Compute autocorrelation
    void Compute(const AllPara&, const OpCorrLocalInfo&);

    // Output results
    void Output(const AllPara&, const OpCorrLocalPara&, const string&);

    ~OpAutoCorr() {};
};

// This will return the parameter value
typedef double (OpCorrLocalPara::*Para_cal)(int, AllPara&);
typedef void (OpCorrLocalPara::*Para_init)(const AllPara&);

// Compute proper local parameters to update in a loop
class OpCorrLocalPara
{
private:
    double min_para_;
    double max_para_;
    int para_pts_; // Total number of changes for parameters

    // The map that records which parameter to update
    map<string, bool> para_update_bool_;

    // The map that calls function for updating parameters
    map<string, Para_cal> para_update_func_;

    // The map that calls function to update min_para_ and max_para_
    map<string, Para_init> para_init_func_;

    string para_name_; // The name of the parameter update that would be called
    string para_; // The name of the parameter

    // Initialize and Update tau_choice
    void Tau_init_(const AllPara&);
    double Tau_update_(int, AllPara&);

public:
    OpCorrLocalPara(const AllPara&);

    double Para_Update(int n, AllPara& local_para) {
        return ( ( this ->*(para_update_func_[para_name_]) )(n, local_para) );
    }
    double Min_Para() const {return min_para_;}
    double Max_Para() const {return max_para_;}
    string Para_Name() const {return para_;}

    ~OpCorrLocalPara() {};
};

#endif //MBL_V1_FLO_OP_AUTO_CORR_H
