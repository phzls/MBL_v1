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
}

