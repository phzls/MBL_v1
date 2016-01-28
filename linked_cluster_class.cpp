//
// Created by Liangsheng Zhang on 1/27/16.
//

#include "linked_cluster_class.h"

using namespace std;

void LinkedCluster::cluster_initialize_(const AllPara& parameters) {

    map<string,bool>::const_iterator it = parameters.linked_cluster_para.linked_cluster_cal.begin();
    while(it != parameters.linked_cluster_para.linked_cluster_cal.end()){
        if(it->second){
            string name = it->first;
            if(name == "Average_Entropy") cluster_cal_.push_back( new ClusterAveEnt(parameters) );
            else{
                cout << "Unknown name " << name << " for linked cluster computation" << endl;
                abort();
            }
        }
    }
}
