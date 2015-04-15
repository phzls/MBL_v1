//
// Created by Liangsheng Zhang on 4/15/15.
//

#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <Eigen/Eigenvalues>
#include "disorder_model_transition.h"
#include "methods.h"
#include "generic_func.h"
#include "screen_output.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains functions related to ent_var, which is the entropy variance for all eigenstates given a
 ** model. The entropy is half chain entanglement entropy, so the spin chain must have even length.
 **/

/*
 * Initialize ent_var in model_data_. The outer index is for J, and the inner index is for realization
 */
void DisorderModelTransition::Ent_var_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.ent_var.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.ent_var[i].resize(num_realization);

        for (int j=0; j<num_realization;j++) model_data_.ent_var[i][j] = 0;
    }

    op_diag_ = true;
    op_dial_evec_keep_ = true;
}

/*
 * Compute ent_var given J index and realization index
 */
void DisorderModelTransition::Ent_var_compute_(AllPara const & parameters, const EvolOP* model,
                                               const DisorderLocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    const int dim = model -> Get_Dim();
    const int size = model -> Get_Size(); // Size of the chain

    vector<double> ent(dim); // Entropy of all eigenstates

    if (size %2 != 0){
        cout << "Length of the chain must be even for entropy variance calculation in flo_transition." << endl;
        cout << "Spin chain length: " << size << endl;
        abort();
    }

    if (local_info.evec_type_real){
        // Vector for eigenvectors in basic binary basis
        vector<vector<double > > evec_basic(dim);
        for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, local_info.evec_real, evec_basic);

        VectorXd evec(dim); // Vector for one eigenstate in binary basis
        MatrixXd reduced_d; // Reduced density matrix

        for (int i=0; i<evec_basic.size(); i++){
            ent[i] = 0;
            for (int j=0; j<evec_basic[i].size();j++) evec[j] = evec_basic[i][j];

            // Compute left reduce density matrix
            reduced_density_left_2(evec, size, size/2, reduced_d);

            SelfAdjointEigenSolver<MatrixXd> density_eigen; // Eigen for reduced density matrix
            density_eigen.compute(reduced_d, false); // Eigenvectors not computed

            // Compute entropy for this state
            for (int j=0; j<density_eigen.eigenvalues().rows();j++){
                double eval = density_eigen.eigenvalues()(j);
                if (abs(eval)>1.0e-10)
                {
                    if (eval<0){
                        cout << "Density matrix has significant negative eigenvalues." << endl;
                        cout << eval << endl;
                        abort();
                    }
                    ent[i] += -eval*log2(eval);
                }
            }

            if (debug){
                cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
                cout << "Reduced Density Matrix: " << endl;
                complex_matrix_write(reduced_d);
                cout << "Entropy: " << ent[i] << endl;
            }
        }
    }
    else{
        // Vector for eigenvectors in basic binary basis
        vector<vector<complex<double> > > evec_basic(dim);
        for (int i=0; i<dim;i++) evec_basic[i].resize(dim);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, local_info.evec_complex, evec_basic);

        VectorXcd evec(dim); // Vector for one eigenstate in binary basis
        MatrixXcd reduced_d; // Reduced density matrix

        for (int i=0; i<evec_basic.size(); i++){
            ent[i] = 0;
            for (int j=0; j<evec_basic[i].size();j++) evec[j] = evec_basic[i][j];

            // Compute left reduce density matrix
            reduced_density_left_2(evec, size, size/2, reduced_d);

            SelfAdjointEigenSolver<MatrixXcd> density_eigen; // Eigen for reduced density matrix
            density_eigen.compute(reduced_d, false); // Eigenvectors not computed

            // Compute entropy for this state
            for (int j=0; j<density_eigen.eigenvalues().rows();j++){
                double eval = density_eigen.eigenvalues()(j);
                if (abs(eval)>1.0e-10)
                {
                    if (eval<0){
                        cout << "Density matrix has significant negative eigenvalues." << endl;
                        cout << eval << endl;
                        abort();
                    }
                    ent[i] += -eval*log2(eval);
                }
            }

            if (debug){
                cout << "Realization " << local_info.realization_index << " eigenstate " << i << ":" << endl;
                cout << "Reduced Density Matrix: " << endl;
                complex_matrix_write(reduced_d);
                cout << "Entropy: " << ent[i] << endl;
            }
        }
    }

    double mean,sd;
    generic_mean_sd(ent, mean, sd);
    model_data_.ent_var[local_info.J_index][local_info.realization_index] = sd;


    if (debug){
        cout << "Realization " << local_info.realization_index << " entropy variance: "
        << sd << endl;
        cout << endl;
    }
}

/*
 * Output averages of mean entropy variance (over all eigenstates for each realization) and their
 * standard deviations among all realizations for each J
 */
void DisorderModelTransition::Ent_var_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.ent_var.size() != J_N){
        cout << "Not enough number of J for entropy variance." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.ent_var.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.ent_var[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for entropy varaince." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.ent_var[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",entropy_variance";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        generic_mean_sd(model_data_.ent_var[i], mean, sd);
        // Convert sd to var
        fout << setw(10) << J << setw(width) << mean << setw(width) << sd * sd << endl;
    }
}




