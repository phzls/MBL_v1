//
// Created by Liangsheng Zhang on 1/27/16.
//

/*
 * This file declares classes used for 1D linked cluster calculations
 */

#ifndef MBL_V1_LINKED_CLUSTER_CLASS_H
#define MBL_V1_LINKED_CLUSTER_CLASS_H

#include <string>
#include <map>
#include <vector>
#include "parameters.h"
#include "evol_op.h"

// A structure containing all information
struct ClusterData
{
    // Model parameters for all realizations
    std::vector< std::vector< std::vector<double> > > model_para;
};

// A structure containing all local information at each step
struct ClusterLocalInfo
{
    int J_index; // The index for J
    int run_index; // The index for realization
};

// A base class for various calculations
class ClusterCal
{
protected:
    string repr_; // Representation string

public:
    virtual void Compute(const AllPara&, const ClusterData&, const ClusterLocalInfo&) = 0;
    // The string gives the first part of filename
    virtual void Output(const AllPara&, const std::string&) = 0;

    // The name of the corresponding cluster computation
    string Name() const {return repr_;}
};



class LinkedCluster
{
private:
    // All possible computations
    std::vector<ClusterCal*> cluster_cal_;

    // Initialize all maps
    void cluster_initialize_(const AllPara&);

public:
    LinkedCluster(const AllPara& parameters){
        cluster_initialize_(parameters);
    }

    void Compute(const AllPara& local_parameters, const ClusterData& cluster_data,
                 const ClusterLocalInfo& local_info){
        for(int i=0; i<cluster_cal_.size(); i++){
            cluster_cal_[i]->Compute(local_parameters, cluster_data, local_info);
        }
    }

    void Output(const AllPara& parameters, const string& name){
        for(int i=0; i<cluster_cal_.size(); i++){
            cluster_cal_[i]->Output(parameters, name);
        }
    }

    ~LinkedCluster(){
        for(int i=0; i<cluster_cal_.size();i++){
            delete cluster_cal_[i];
            cluster_cal_[i] = NULL;
        }
    };
};

// Average entropy computation
class ClusterAveEnt : public ClusterCal
{
private:
    // A vector to record entropy. The inner index is for different realization.
    // The middle index is for different orders
    // The outer index is for different J
    std::vector< std::vector< std::vector<double> > > ave_ent_;

    const int order_N; // Number of orders
    const int J_N; // Number of J
    const int run_N; // Number of realizations
    const string model_name; // Name of the model
    const bool debug; // Whether the program will output debug information
public:
    // std::vector<EvolOP*> op_vec; // A vector of various operators
    void Compute(const AllPara&, const ClusterData&, const ClusterLocalInfo&);
    void Output(const AllPara&, const string&);

    ClusterAveEnt(const AllPara&);
    ~ClusterAveEnt(){};
};

#endif //MBL_V1_LINKED_CLUSTER_CLASS_H
