//
// Created by Liangsheng Zhang on 1/27/16.
//

/*
 * This file computes average entropy across a fixed cut for 1D linked clusters
 * Here the averaging refers to averaging over different eigenstates. The cut
 * is put at the middle
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <iomanip>
#include <Eigen/Eigenvalues>
#include "linked_cluster_class.h"
#include "methods.h"
#include "generic_func.h"
#include "tasks_models.h"

using namespace std;
using namespace Eigen;

extern TasksModels tasks_models; // Record all the tasks and methods. Defined in main.

// Initialize the class
ClusterAveEnt::ClusterAveEnt(const AllPara& parameters) :
        order_N(parameters.linked_cluster_para.order),
        J_N(parameters.floquet.J_N), run_N(parameters.generic.num_realizations),
        model_name(parameters.generic.model), debug(parameters.generic.debug)
{
    if(order_N<2){
        cout << "For linked cluster calculation of entropy order must be larger than 2" << endl;
        cout << "Order: "<< order_N << endl;
        abort();
    }

    ave_ent_.resize(J_N);
    for(int i=0; i<J_N; i++){
        ave_ent_[i].resize(order_N-1);

        for(int j=0; j<order_N-1; j++){
            ave_ent_[i][j].resize(run_N);

            for(int k=0; k<run_N; k++) ave_ent_[i][j][k] = 0;
        }
    }

    repr_ = "Average Entropy";
    model_type_determined_ = false;
}

// Compute entropies from all orders by first constructing models.
void ClusterAveEnt::Compute(const AllPara& parameters, const ClusterData& cluster_data,
                            const ClusterLocalInfo& local_info) {
    // Total number of sites used on the maximum chain. The cut is at the middle
    const int total_size = 2*(order_N-1);

    for(int i=0; i<cluster_data.model_para[local_info.run_index].size(); i++){
        if(cluster_data.model_para[local_info.run_index][i].size() != total_size){
            cout << "Incompatible order to number of model parameters at " << local_info.run_index <<
                    "th realization and " << i << "th vector" << endl;
            cout << "Number of order: " << order_N << endl;
            cout << "Number of model parameters: " <<
                    cluster_data.model_para[local_info.run_index][i].size() << endl;
            abort();
        }
    }

    AllPara local_parameters(parameters); // Local parameters which can be changed

    // Record of all entropy. The outer index is order, and the inner index is the starting
    // position of the sub-chain
    vector< vector<double> > all_ent(order_N-1);

    if(debug) cout << "J index: " << local_info.J_index << " run index: "
                   << local_info.run_index << endl;
    for(int i=0; i<order_N-1; i++){
        int order = i+2; // This defines the length of the model
        local_parameters.generic.size = order; // Size of the system
        local_parameters.floquet.J = local_info.J; // Coupling strength

        all_ent[i].resize(order-1);

        // Model parameters
        vector< vector<double> > model_para(cluster_data.model_para[local_info.run_index].size(),
                                            vector<double>(order));

        // For starting position on the maximum chain
        for(int j= order_N-order; j<order_N-1; j++){
            // Initialize parameters
            for(int k=0; k<model_para.size(); k++){
                for(int l=0; l<order; l++){
                    model_para[k][l] = cluster_data.model_para[local_info.run_index][k][j+l];
                }
            }

            if(debug) cout << "order: " << i <<  " starting position: " << j << endl;

            EvolOP* model;
            tasks_models.Model(model_name, local_parameters, model);
            model->Evol_Para_Copy(model_para); // Copy parameters
            model->Evol_Construct(); // Construct the model matrix
            model->Evol_Diag(true); // Diagonalize the matrix and keep eigenvectors
            model->OP_Erase(); // Erase the matrix

            if(local_info.J_index == 0 && local_info.run_index == 0 && !model_type_determined_){
                model_type_ = model->Type();
                model_type_determined_ = true;
            }

            vector<double> ent_vec; // Record entropy for each eigenstate

            model_entropy_left_2(model, order_N-1-j, ent_vec); // Compute entropy for all eigenstates

            delete model;
            model = NULL;

            double ave_ent = 0; // Average entropy for all eigenstates
            for(int l=0; l<ent_vec.size();l++) ave_ent += ent_vec[l];
            ave_ent /= double(ent_vec.size());

            if(debug){
                cout << "Raw entropy: " << ave_ent << endl;
            }


            // Subtract entropy of all embedied sub-chains
            for(int l=0; l<i; l++){ // The order of the sub-chain
                int l_order = l+2;
                int min = j; // Minimum starting position
                if(min < order_N-l_order) min = order_N - l_order;

                int max = j+order-l_order; // Maximum starting position
                if(max>order_N-2) max = order_N-2;

                for(int k=min; k<=max; k++){
                    ave_ent -= all_ent[l][k-order_N+l_order];
                }
            }

            if(debug) cout << "Linked cluster entropy: " << ave_ent << endl;

            all_ent[i][j-order_N+order] = ave_ent;
        }

        ave_ent_[local_info.J_index][i][local_info.run_index] = 0;
        for(int j=0; j<all_ent[i].size(); j++){
            ave_ent_[local_info.J_index][i][local_info.run_index] += all_ent[i][j];
        }

        if(debug) cout << "order: " <<i << " ave ent: "
                       << ave_ent_[local_info.J_index][i][local_info.run_index] << endl;
    }
}

// Output the average and unbiased standard deviation of coefficients of all orders
// from all realizations
void ClusterAveEnt::Output(const AllPara& parameters, const string& name) {
    const double J_min = parameters.floquet.J_min;
    const double J_max = parameters.floquet.J_max;
    const int version = parameters.generic.version;
    const int width = parameters.output.width;
    const bool output = parameters.output.filename_output;
    const int size = parameters.generic.size;

    if (ave_ent_.size() != J_N){
        cout << "Not enough number of J for linked cluster average entropy." << endl;
        cout << "Expected Number: " << J_N << endl;
        cout << "Actual number: " << ave_ent_.size() << endl;
        abort();
    }

    for (int i=0; i< J_N; i++){
        if (ave_ent_[i].size() != order_N-1){
            cout << "Wrong number of order at " << i <<"th J for linked cluster average entropy."
                 << endl;
            cout << "Expected Number: " << order_N-1 << endl;
            cout << "Actual Number: " << ave_ent_[i].size() << endl;
            abort();
        }
        for(int j=0; j<order_N-1; j++){
            if (ave_ent_[i][j].size() != run_N){
                cout << "Not enough number of realizations at " << i <<"th J "
                     << "and " << j << "th order for linked cluster average entropy." << endl;
                cout << "Expected Number: " << run_N << endl;
                cout << "Actual Number: " << ave_ent_[i][j].size() << endl;
                abort();
            }
        }
    }

    stringstream filename;
    string new_name = name;
    if(name.size() > 0) new_name += ",";
    filename << new_name << model_type_ << ",size=" << size << ",Run=" << run_N
    << ",J_N=" << J_N << ",J_min=" << J_min << ",J_max=" << J_max <<",order="
    << order_N << ",ave_entropy";
    if (version > 0) filename <<",v" << version;
    filename << ".txt";

    if (output) cout << filename.str() <<endl;

    ofstream fout( filename.str().c_str() );

    for (int i=0; i<J_N;i++){
        double J;
        if (J_N > 1) J = J_min + i * (J_max - J_min)/(J_N-1);
        else J = J_min;

        for(int j=0; j<order_N-1;j++){
            double order = j+2;
            double mean,sd;
            generic_mean_sd(ave_ent_[i][j], mean, sd);

            // Convert sd to var
            fout << setw(10) << J << setw(width) << order << setw(width) << mean
                 << setw(width) << sd << endl;
        }
    }
}
