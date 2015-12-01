//
// Created by Liangsheng Zhang on 4/18/15.
//

/*
 * Functions related to reading parameters from files for all models and tasks
 */
#include "para_model_task.h"

using namespace std;

//=================================== MODELS ===================================================

// For Ising random simple floquet operator
void Ising_Random_Simp_Flo_Para(AllPara& parameters, string count){
    string filename = "ising_random_simp_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Ising random simple shift real floquet operator
void Ising_Random_Simp_Shift_Real_Flo_Para(AllPara& parameters, string count){
    string filename = "ising_random_simp_shift_real_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Ising quasi-periodic simple shift real floquet operator
void Ising_Quasi_Simp_Shift_Real_Flo_Para(AllPara& parameters, string count){
    string filename = "ising_quasi_simp_shift_real_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Ising all random simple shift real floquet operator
void Ising_All_Random_Simp_Shift_Real_Flo_Para(AllPara& parameters, string count){
    string filename = "ising_all_random_simp_shift_real_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Ising all quasi-periodic simple shift real floquet operator
void Ising_All_Quasi_Simp_Shift_Real_Flo_Para(AllPara& parameters, string count){
    string filename = "ising_all_quasi_simp_shift_real_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Ising random simple shift cosine real floquet operator
void Ising_Random_Simp_Shift_Cos_Real_Flo_Para(AllPara& parameters, string count){
    string filename = "ising_random_simp_shift_cos_real_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Heisenberg random cosine Sz sector Hamiltonian Operator
void Heisen_Random_Cos_Sz_Sector_Ham_Para(AllPara& parameters, string count){
    string filename = "heisen_random_cos_sz_sector_ham_" + count;
    vector<vector<string> > content;
    string keyword = "J";
    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Heisenberg quasi-periodic cosine Sz sector Hamiltonian Operator
void Heisen_Quasi_Sz_Sector_Ham_Para(AllPara& parameters, string count){
    string filename = "heisen_quasi_sz_sector_ham_" + count;
    vector<vector<string> > content;
    string keyword = "J";
    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Ising random simple shift cosine real tau floquet operator
void Ising_Random_Simp_Shift_Cos_Real_Tau_Flo_Para(AllPara& parameters, string count){
    string filename = "ising_random_simp_shift_cos_real_tau_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "tau";
    para_get(filename, content, keyword, parameters.floquet.tau);
}

// For Heisenberg random cos sz sector shift real tau floquet operator
void Heisen_Random_Cos_Sz_Sector_Shift_Real_Tau_Flo_Para(AllPara& parameters, string count){
    string filename = "heisen_random_cos_sz_sector_shift_real_tau_flo_" + count;
    vector< vector<string> > content;
    para_file_read(filename, content);

    string keyword = "J";
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "tau";
    para_get(filename, content, keyword, parameters.floquet.tau);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Heisenberg quasi-periodic sz sector shift real tau floquet operator
void Heisen_Quasi_Sz_Sector_Shift_Real_Tau_Flo_Para(AllPara& parameters, string count){
    string filename = "heisen_quasi_sz_sector_shift_real_tau_flo_" + count;
    vector< vector<string> > content;
    para_file_read(filename, content);

    string keyword = "J";
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "tau";
    para_get(filename, content, keyword, parameters.floquet.tau);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Ising random simple cosine hamiltonian operator
void Ising_Random_Simp_Cos_Ham_Para(AllPara& parameters, string count){
    string filename = "ising_random_simp_cos_ham_" + count;
    vector< vector<string> > content;
    para_file_read(filename, content);

    string keyword = "J";
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Ising quasi-periodic simple hamiltonian operator
void Ising_Quasi_Simp_Ham_Para(AllPara& parameters, string count){
    string filename = "ising_quasi_simp_ham_" + count;
    vector< vector<string> > content;
    para_file_read(filename, content);

    string keyword = "J";
    para_get(filename, content, keyword, parameters.floquet.J);
}

// For Heisenberg random cos sz sector modified tau floquet operator
void Heisen_Random_Cos_Sz_Sector_Modified_Tau_Flo_Para(AllPara& parameters, string count){
    string filename = "heisen_random_cos_sz_sector_modified_tau_flo_" + count;
    vector< vector<string> > content;
    para_file_read(filename, content);

    string keyword = "J";
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "tau";
    para_get(filename, content, keyword, parameters.floquet.tau);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Heisenberg quasi-periodic sz sector modified tau floquet operator
void Heisen_Quasi_Sz_Sector_Modified_Tau_Flo_Para(AllPara& parameters, string count){
    string filename = "heisen_quasi_sz_sector_modified_tau_flo_" + count;
    vector< vector<string> > content;
    para_file_read(filename, content);

    string keyword = "J";
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "tau";
    para_get(filename, content, keyword, parameters.floquet.tau);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Heisenberg modified random cosine Sz sector Hamiltonian Operator
void Heisen_Modified_Random_Cos_Sz_Sector_Ham_Para(AllPara& parameters, string count){
    string filename = "heisen_modified_random_cos_sz_sector_ham_" + count;
    vector<vector<string> > content;
    string keyword = "J";
    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Heisenberg modified quasi-periodic cosine Sz sector Hamiltonian Operator
void Heisen_Modified_Quasi_Sz_Sector_Ham_Para(AllPara& parameters, string count){
    string filename = "heisen_modified_quasi_sz_sector_ham_" + count;
    vector<vector<string> > content;
    string keyword = "J";
    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);
}

// For Heisenberg continuous modified random cosine Sz sector Hamiltonian Operator
void Heisen_Con_Modified_Random_Cos_Sz_Sector_Ham_Para(AllPara& parameters, string count){
    string filename = "heisen_con_modified_random_cos_sz_sector_ham_" + count;
    vector<vector<string> > content;
    string keyword = "J";
    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);

    keyword = "alpha";
    para_get(filename, content, keyword, parameters.floquet.alpha);
}

// For Heisenberg continuous modified quasi-periodic cosine Sz sector Hamiltonian Operator
void Heisen_Con_Modified_Quasi_Sz_Sector_Ham_Para(AllPara& parameters, string count){
    string filename = "heisen_con_modified_quasi_sz_sector_ham_" + count;
    vector<vector<string> > content;
    string keyword = "J";
    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);

    keyword = "Total_Spin_Z";
    para_get(filename, content, keyword, parameters.floquet.total_spin_z);

    keyword = "alpha";
    para_get(filename, content, keyword, parameters.floquet.alpha);
}

// For XXZ Gaussian random shift real floquet operator
void XXZ_Gaussian_Random_Shift_Real_Flo_Para(AllPara& parameters, string count){
    string filename = "xxz_gaussian_random_shift_real_flo_" + count;
    vector<vector<string> > content;
    string keyword = "J";

    para_file_read(filename, content);
    para_get(filename, content, keyword, parameters.floquet.J);
}



//=================================== TASKS =====================================================

// generic parameters
void Generic_Para(AllPara& parameters, string count){
    string filename = "generic_para_" + count;
    vector<vector<string> > content;

    para_file_read(filename, content);

    string keyword = "task";
    para_get(filename, content, keyword, parameters.generic.task);

    keyword = "model";
    para_get(filename, content, keyword, parameters.generic.model);

    keyword = "size";
    para_get(filename, content, keyword, parameters.generic.size);

    keyword = "num_realizations";
    para_get(filename, content, keyword, parameters.generic.num_realizations);

    keyword = "threads_N";
    para_get(filename, content, keyword, parameters.generic.threads_N);

    keyword = "debug";
    para_get(filename, content, keyword, parameters.generic.debug);

    keyword = "version";
    para_get(filename, content, keyword, parameters.generic.version);
}

// For output parameters
void Output_Para(AllPara& parameters, string count){
    string filename = "output_para_" + count;
    vector<vector<string> > content;

    para_file_read(filename, content);

    string keyword = "time";
    para_get(filename, content, keyword, parameters.generic.time);

    keyword = "width";
    para_get(filename, content, keyword, parameters.output.width);

    keyword = "filename_output";
    para_get(filename, content, keyword, parameters.output.filename_output);
}

// For disorder_transition
void Disorder_Transition_Para(AllPara& parameters, string count){
    string filename = "disorder_transition_" + count;
    vector<vector<string> > content;

    para_file_read(filename, content);

    string keyword = "J_N";
    para_get(filename, content, keyword, parameters.floquet.J_N);

    keyword = "J_max";
    para_get(filename, content, keyword, parameters.floquet.J_max);

    keyword = "J_min";
    para_get(filename, content, keyword, parameters.floquet.J_min);

    keyword = "mid_half_spectrum";
    para_get(filename, content, keyword, parameters.transition.mid_half_spectrum);

    keyword = "ZZ_Correlation_Square";
    bool choice;
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["ZZ_Correlation_Square"] = choice;

    keyword = "Entropy_Variance";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["Entropy_Variance"] = choice;

    keyword = "ZZ_Time_Correlation";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["ZZ_Time_Correlation"] = choice;

    keyword = "ZZ_Time_Correlation_Components";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["ZZ_Time_Correlation_Components"] = choice;

    keyword = "ZZ_All_Correlation_Square";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["ZZ_All_Correlation_Square"] = choice;

    keyword = "ZZ_Correlation_Square_All_Sample";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["ZZ_Correlation_Square_All_Sample"] = choice;

    keyword = "ZZ_Time_Correlation_All_Sample";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["ZZ_Time_Correlation_All_Sample"] = choice;

    keyword = "Log_ZZ_Correlation_Square";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["Log_ZZ_Correlation_Square"] = choice;

    keyword = "Log_ZZ_All_Correlation_Square";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["Log_ZZ_All_Correlation_Square"] = choice;

    keyword = "ZZ_All_Time_Correlation";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["ZZ_All_Time_Correlation"] = choice;

    keyword = "Entropy_Variance_All_Mean";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute["Entropy_Variance_All_Mean"] = choice;

    keyword = "Flo_Level_Stats_Ave";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute[keyword] = choice;

    keyword = "Ent_Scaled_Mean";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute[keyword] = choice;

    keyword = "Ham_Level_Stats_Ave";
    para_get(filename, content, keyword, choice);
    parameters.transition.flo_transition_compute[keyword] = choice;
}

// For single_model_time_evolution_para
void Single_Model_Time_Evolution_Para(AllPara& parameters, string count){
    string filename = "single_model_time_evolution_" + count;
    vector<vector<string> > content;

    para_file_read(filename, content);

    string keyword = "time_step";
    para_get(filename, content, keyword, parameters.evolution.time_step);

    keyword = "step_size";
    para_get(filename, content, keyword, parameters.evolution.step_size);

    keyword = "init_func_name";
    para_get(filename, content, keyword, parameters.evolution.init_func_name);

    keyword = "log_time";
    para_get(filename, content, keyword, parameters.evolution.log_time);

    keyword = "log_time_jump";
    para_get(filename, content, keyword, parameters.evolution.log_time_jump);

    keyword = "jump";
    para_get(filename, content, keyword, parameters.evolution.jump);

    keyword = "Entropy_Per_Model";
    bool choice;
    para_get(filename, content, keyword, choice);
    parameters.evolution.evol_compute["Entropy_Per_Model"] = choice;

    keyword = "Leftmost_Spin_Per_Model";
    choice;
    para_get(filename, content, keyword, choice);
    parameters.evolution.evol_compute["Leftmost_Spin_Per_Model"] = choice;

    keyword = "left_size";
    para_get(filename, content, keyword, parameters.evolution.left_size);

    keyword = "sample_detail";
    para_get(filename, content, keyword, parameters.evolution.sample_detail);

    keyword = "evol_way";
    para_get(filename, content, keyword, parameters.evolution.evol_way);

    parameters.evolution.model_num = 1;
    parameters.evolution.markov_jump = false;
}

// For multi_model_time_evolution_para
void Multi_Model_Time_Evolution_Para(AllPara& parameters, string count){
    string filename = "multi_model_time_evolution_" + count;
    vector<vector<string> > content;

    para_file_read(filename, content);

    string keyword = "time_step";
    para_get(filename, content, keyword, parameters.evolution.time_step);

    keyword = "step_size";
    para_get(filename, content, keyword, parameters.evolution.step_size);

    keyword = "init_func_name";
    para_get(filename, content, keyword, parameters.evolution.init_func_name);

    keyword = "log_time";
    para_get(filename, content, keyword, parameters.evolution.log_time);

    keyword = "log_time_jump";
    para_get(filename, content, keyword, parameters.evolution.log_time_jump);

    keyword = "jump";
    para_get(filename, content, keyword, parameters.evolution.jump);

    keyword = "Entropy_Per_Model";
    bool choice;
    para_get(filename, content, keyword, choice);
    parameters.evolution.evol_compute["Entropy_Per_Model"] = choice;

    keyword = "Leftmost_Spin_Per_Model";
    choice;
    para_get(filename, content, keyword, choice);
    parameters.evolution.evol_compute["Leftmost_Spin_Per_Model"] = choice;

    keyword = "left_size";
    para_get(filename, content, keyword, parameters.evolution.left_size);

    keyword = "sample_detail";
    para_get(filename, content, keyword, parameters.evolution.sample_detail);

    keyword = "model_num";
    para_get(filename, content, keyword, parameters.evolution.model_num);

    keyword = "evol_way";
    para_get(filename, content, keyword, parameters.evolution.evol_way);

    parameters.evolution.markov_jump = false;
}

// For operator autocorrelation
void Op_Auto_Corr_Para(AllPara& parameters, string count){
    string filename = "op_auto_corr_" + count;
    vector<vector<string> > content;

    para_file_read(filename, content);

    string keyword = "time_pts";
    para_get(filename, content, keyword, parameters.flo_op_auto_corr.time_pts);

    keyword = "para_pts";
    para_get(filename, content, keyword, parameters.flo_op_auto_corr.para_pts);

    keyword = "tau_choice";
    para_get(filename, content, keyword, parameters.flo_op_auto_corr.tau_choice);

    keyword = "tau_max";
    para_get(filename, content, keyword, parameters.flo_op_auto_corr.tau_max);

    keyword = "tau_min";
    para_get(filename, content, keyword, parameters.flo_op_auto_corr.tau_min);

    vector<string> op_corr_map;
    vector<string>::iterator it;

    op_corr_map.push_back("Energy");

    for(it = op_corr_map.begin(); it != op_corr_map.end(); it++){
        bool choice;
        keyword = *it;
        para_get(filename, content, keyword, choice);
        parameters.flo_op_auto_corr.op_corr_map[keyword] = choice;
    }
}

