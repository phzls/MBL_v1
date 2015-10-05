//
// Created by Liangsheng Zhang on 9/16/15.
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
#include "ham_evol_model.h"
#include "flo_evol_model.h"

using namespace std;
using namespace Eigen;

/**
 ** This file contains function related to ent_scaled_mean, which the average of entanglement entropy
 ** for all (chosen) eigenstates for each sample, scaled by Page value up to first correction to L/2.
 **/

/*
 * Initialize ent_scaled_mean in model_data_. The outer index is for J, and the inner index is for realization
 */
void DisorderModelTransition::Ent_scaled_mean_init_(AllPara const & parameters) {
    const int J_N = parameters.floquet.J_N; // Number of points for J
    const int num_realization = parameters.generic.num_realizations;

    model_data_.ent_scaled_mean.resize(J_N);

    for (int i=0; i<J_N; i++){
        model_data_.ent_scaled_mean[i].resize(num_realization);

        for (int j=0; j<num_realization;j++) model_data_.ent_scaled_mean[i][j] = 0;
    }

    op_diag_ = true;
    op_dial_evec_keep_ = true;
}

/*
 * Compute ent_scaled_mean given J index and realization index
 */
void DisorderModelTransition::Ent_scaled_mean_compute_(AllPara const & parameters, const EvolOP* model,
                                                       const DisorderLocalInfo& local_info) {

    const bool debug = parameters.generic.debug;

    const int dim = model -> Get_Dim();
    const int size = model -> Get_Size(); // Size of the chain

    const double Page_val = size/2 - 1.0/2.0/log(2);

    vector<double> ent(dim); // Entropy of all eigenstates

    if (size %2 != 0){
        cout << "Length of the chain must be even for entropy variance calculation in flo_transition." << endl;
        cout << "Spin chain length: " << size << endl;
        abort();
    }

    if ( model -> Eigen_Basis_Type() == "Restricted Basic" && local_info.evec_type_real ){
        // This is just a temporary workaround

        int half_mid_size = dim - 2*(dim/4);

        // May choose to use all eigenstates
        if(!local_info.is_ham || !local_info.mid_half_spectrum) half_mid_size = dim;

        ent.resize(half_mid_size);

        int full_dim = 1 << size; // Full dimension of the spin chain

        vector<int> spin_config;
        VectorXd evec(full_dim); // Used to store one eigenvector
        MatrixXd reduced_d; // Reduced density matrix

        if( dynamic_cast<const HamEvolHeisenRandomCosSzSector*>(model) ){
            const HamEvolHeisenRandomCosSzSector* down_cast_model = dynamic_cast<const HamEvolHeisenRandomCosSzSector*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }
        else if( dynamic_cast<const HamEvolHeisenQuasiSzSector*>(model) ){
            const HamEvolHeisenQuasiSzSector* down_cast_model = dynamic_cast<const HamEvolHeisenQuasiSzSector*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }
        else if(dynamic_cast<const FloEvolHeisenRandomCosSzSectorShiftRealTau*>(model) ){
            const FloEvolHeisenRandomCosSzSectorShiftRealTau* down_cast_model =
                    dynamic_cast<const FloEvolHeisenRandomCosSzSectorShiftRealTau*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }
        else if(dynamic_cast<const FloEvolHeisenQuasiSzSectorShiftRealTau*>(model) ){
            const FloEvolHeisenQuasiSzSectorShiftRealTau* down_cast_model =
                    dynamic_cast<const FloEvolHeisenQuasiSzSectorShiftRealTau*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }
        else if(dynamic_cast<const FloEvolHeisenRandomCosSzSectorModifiedTau*>(model) ){
            const FloEvolHeisenRandomCosSzSectorModifiedTau* down_cast_model =
                    dynamic_cast<const FloEvolHeisenRandomCosSzSectorModifiedTau*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }
        else if(dynamic_cast<const FloEvolHeisenQuasiSzSectorModifiedTau*>(model) ){
            const FloEvolHeisenQuasiSzSectorModifiedTau* down_cast_model =
                    dynamic_cast<const FloEvolHeisenQuasiSzSectorModifiedTau*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }
        else if( dynamic_cast<const HamEvolHeisenModifiedRandomCosSzSector*>(model) ){
            const HamEvolHeisenModifiedRandomCosSzSector* down_cast_model =
                    dynamic_cast<const HamEvolHeisenModifiedRandomCosSzSector*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }
        else if( dynamic_cast<const HamEvolHeisenModifiedQuasiSzSector*>(model) ){
            const HamEvolHeisenModifiedQuasiSzSector* down_cast_model =
                    dynamic_cast<const HamEvolHeisenModifiedQuasiSzSector*>(model);
            spin_config = down_cast_model -> Get_Spin_Config();
        }

        int index = 0;
        for(int i=0; i<local_info.evec_real.size(); i++){
            for (int j=0; j<local_info.evec_real[i].cols(); j++){
                if(index >= half_mid_size){
                    cout << "Too many eigenvectors for middle half spectrum." << endl;
                    cout << "Expected: " << half_mid_size << endl;
                    cout << "Current index: " << index << endl;
                    abort();
                }

                ent[index] = 0;
                evec.setZero();

                for(int k=0; k<local_info.evec_real[i].rows(); k++){
                    int pos = spin_config[k];
                    evec(pos) = local_info.evec_real[i](k,j);
                }

                // Compute left reduce density matrix
                reduced_density_left_2(evec, size, size/2, reduced_d);

                SelfAdjointEigenSolver<MatrixXd> density_eigen; // Eigen for reduced density matrix
                density_eigen.compute(reduced_d, false); // Eigenvectors not computed

                // Compute entropy for this state
                for (int l=0; l<density_eigen.eigenvalues().rows();l++){
                    double eval = density_eigen.eigenvalues()(l);
                    if (abs(eval)>1.0e-10)
                    {
                        if (eval<0){
                            cout << "Density matrix has significant negative eigenvalues." << endl;
                            cout << eval << endl;
                            abort();
                        }
                        ent[index] += -eval*log2(eval);
                    }
                }

                if (debug){
                    cout << "Realization " << local_info.realization_index << " eigenstate " << index << ":" << endl;
                    cout << "Full Eigenvector:" << endl;
                    cout << endl;
                    real_matrix_write(evec);
                    cout << "Reduced Density Matrix: " << endl;
                    real_matrix_write(reduced_d);
                    cout << "Entropy: " << ent[index] << endl;
                }
                ++ index;
            }
        }

        if(index != half_mid_size){
            cout << "Incorrect number of eigenstates in entropy calculation." << endl;
            cout << "Expected: " << half_mid_size << endl;
            cout << "Obtained: " << index << endl;
            abort();
        }
    }
    else if (local_info.evec_type_real){
        // Vector for eigenvectors in basic binary basis
        vector<vector<double > > evec_basic(dim);

        // For half spectrum in Hamiltonian
        if(local_info.is_ham && local_info.mid_half_spectrum){
            int half_mid_size = dim - 2*(dim/4);
            evec_basic.resize(half_mid_size);
            ent.resize(half_mid_size);
        }
        for (int i=0; i<evec_basic.size();i++) evec_basic[i].resize(dim);

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
                real_matrix_write(reduced_d);
                cout << "Entropy: " << ent[i] << endl;
            }
        }
    }
    else{
        // Vector for eigenvectors in basic binary basis
        vector<vector<complex<double> > > evec_basic(dim);

        // For half spectrum in Hamiltonian
        if(local_info.is_ham && local_info.mid_half_spectrum){
            int half_mid_size = dim - 2*(dim/4);
            evec_basic.resize(half_mid_size);
            ent.resize(half_mid_size);
        }
        for (int i=0; i<evec_basic.size();i++) evec_basic[i].resize(dim);

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
    model_data_.ent_scaled_mean[local_info.J_index][local_info.realization_index] = mean / Page_val;

    if (debug){
        cout << "Realization " << local_info.realization_index << " scaled average entropy: "
             << model_data_.ent_scaled_mean[local_info.J_index][local_info.realization_index] << endl;
    }
}

/*
 * Output averages of scaled mean entropy (over all eigenstates for each realization) and their
 * standard deviations among all realizations for each J
 */
void DisorderModelTransition::Ent_scaled_mean_out_(AllPara const & parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int J_N = parameters.floquet.J_N;
    const int num_realizations = parameters.generic.num_realizations;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (model_data_.ent_scaled_mean.size() != J_N){
        cout << "Not enough number of J for scaled entropy mean." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << model_data_.ent_scaled_mean.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (model_data_.ent_scaled_mean[i].size() != num_realizations){
            cout << "Not enough number of realizations at " << i <<"th J for scaled entropy mean." << endl;
            cout << "Expected Number: " << num_realizations << endl;
            cout << "Actual Number: " << model_data_.ent_scaled_mean[i].size() << endl;
            abort();
        }
    }

    stringstream filename;
    filename << name << ",size=" << size << ",Run=" << num_realizations << ",J_N=" << J_N << ",J_min=" << J_min
    << ",J_max=" << J_max << ",scaled_entropy_mean";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        double mean,sd;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        generic_mean_sd(model_data_.ent_scaled_mean[i], mean, sd);
        fout << setw(10) << J << setw(width) << mean << setw(width) << sd << endl;
    }
}