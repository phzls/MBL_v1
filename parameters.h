//
// Created by Liangsheng Zhang on 4/15/15.
//

#ifndef MBL_V1_PARAMETERS_H
#define MBL_V1_PARAMETERS_H

/*
 * This file constains structs that give all necessary parameters to pass to various
 * functions. They may contain more parameters than needed for generality.
 */

#include <string>
#include <vector>
#include <utility>
#include <map>

using namespace std;

/*
 * Generic parameters which may be used by any method/function
 */
struct GenericPara
{
	int num_realizations; // Number of realizations
	int threads_N; // Number of threads
	int size; // System size
	int local_dimension; // Local dimension at each site
	bool debug; // Whether show debug information
	bool iso_keep; // Whether isolated part is kept
	string task; // String for the computation task
	string model; // String for the model
	int version; // Version number of the output file. If it is smaller than 1,
				 // then no version number is outputted
	bool time; // Whether various parts of the program is timed
};

/*
 * Parameters that are related to output and writing to files
 */
struct OutputPara
{
	int width; // Width of output in files
	bool filename_output; // Whether write out filenames
};

/*
 * Parameters that are related to Floquet models.
 * For now they are also used for Hamiltonian models
 */
struct FloPara
{
	double J_min; // Minimum J in a loop
	double J_max; // Maximum J in a loop
	int J_N; // Number points for different coupling strength J, including J_min
			 // and J_max
	double J; // Coupling strength J
	double tau; // Time step

	int total_spin_z; // Total Z spin for sector

	double alpha; // The parameter used in continuous modified Hamiltonian
};

struct MatrixPara
{
	string type; // Determine the representation type of the matrix
};

/*
 * Parameters used in time evolutions
 */
struct Evolution
{
	int time_step; // Number of time steps
	double step_size; // Size of time step; would be tau for floquet models
	int model_num; // Number of different models used for time evolution
	int jump; // jump of time points in evolution
	string init_func_name; // Initial state construction function name
	map<string,bool> evol_compute; // Determine which data to compute
	map<string,bool> evol_total_compute; // Determine which data to compute which contains
										 // results for all models
	int left_size; // If partition the chain to two halves, the size of left part

	bool log_time; // Determine whether time increases logarithmically
	int log_time_jump; // The base number which time increases on when logarithmically

	bool markov_jump; // Determine whether there will be markov_time_jump
	int markov_time_jump; // The jump time in markov time evolution

	int leftmost_spin_z_index; // The number gives the index of leftmost spin z value

	bool sample_detail; // Whether output sample_to_sample detail data. If only one model is computed, sample
                        // refers to different initial condition; otherwise it refers to different realization
                        // of model

	string evol_way; // How evolution is done. For now it can only be "vector" or "matrix"
};

/*
 * Parameters used for Markov models
 */
struct Markov
{
	double K; // Coupling strength to the bath
};

/*
 * Initial information related to multiple sets of parameters under one model of evolution
 * Since so far each initial state construction function only requires one set of parameters,
 * this vector should not be used in any of these functions. Instead, the corresponding one-
 * set-parameter variable should be passed in.
 */
struct MultipleInitPara
{
	// The set of numbers which give indices of leftmost spin z value
	vector<int> leftmost_spin_z_index_set;

	// Whether output full evolution results or just some of them from leftmost_spin_z_index_set.
	// If it is true, then only some are output, and whether they change signs are outputted
	bool easy_full_leftmost_spin_z;

	// A small value beyond which is considered as non-zero
	double non_zero_threshold;
};

/*
 * Parameters used for single model computation
 */
struct SingleModel
{
	map<string,bool> single_model_compute; // Determine which method to be computed
};

/*
 * Parameters used for studying transition of models
 */
struct Transition
{
	bool mid_half_spectrum; // Whether only using middle half spectrum for Hamiltonian systems
	map<string,bool> flo_transition_compute; // Determine which method to be computed for floquet systems
};

/*
 * Parameters used for studying eigenstates related quantities
 */
struct Eigenvec
{
	map<string, bool> flo_eigen_compute; // Determine which method to be computed for floquet systems
};

/*
 * Parameters used for studying operator autocorrelation under Floquet dynamics
 */
struct FloOpAutoCorr
{
    bool tau_choice; // Whether tau is varied in studying autocorrelation
    double tau_min;
    double tau_max;
    int para_pts; // Number of points for the changing variable
    int time_pts; // Number of time points

    map<string, bool> op_corr_map; // Determines which operators are used for autocorrelation calculation
};

/*
 * Parameters used for linked cluster calculation
 */
struct LinkedClusterPara
{
    map<string, bool> linked_cluster_cal; // Determine what to compute for linked clusters
    int order;
};

/*
 * All parameters.
 */
struct AllPara
{
	// Generic parameters
	GenericPara generic;

	// Output parameters
	OutputPara output;

	// Floquet parameters
	FloPara floquet;

	// Matrix relevant parameters
	MatrixPara matrix_para;

	// Time evolution parameters
	Evolution evolution;

	// Markov model parameters
	Markov markov;

	// Multiple sets of initial conditions
	MultipleInitPara multi_ini_para;

	// Single model computation
	SingleModel single_model;

	// Model transition computation
	Transition transition;

	// Eigenstates related properties
	Eigenvec eigenvec;

    // Parameters for operator autocorrelation under Floquet dynamics
    FloOpAutoCorr flo_op_auto_corr;

    // Parameters for linked cluster calculation
    LinkedClusterPara linked_cluster_para;
};

#endif //MBL_V1_PARAMETERS_H
