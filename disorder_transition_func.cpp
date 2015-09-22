//
// Created by Liangsheng Zhang on 9/18/15.
//

#include <iostream>
#include <vector>
#include "disorder_transition_func.h"
#include "screen_output.h"

using namespace std;

void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::VectorXd& vec){
    const int size = pair_index.size();
    const int mid_half_size = size - 2*(size/4); // It is at least size/2
    std::vector<bool> index_bool(size, false);

    // Sort the vector according to pair_index
    for(int index=0; index<size;){
        int pos = pair_index[index].second;

        swap(vec(index), vec(pos));
        index_bool[index] = true;

        int next = pair_index[pos].second;

        while(next != index){
            swap(vec(pos), vec(next));
            index_bool[pos] = true;

            pos = next;
            next = pair_index[pos].second;
        }

        index_bool[pos] = true;
        while(index < size && index_bool[index]) ++index;
    }

    // Move middle half to first part
    for(int index=0; index<mid_half_size;index++){
        swap(vec(index), vec(index+size/4));
    }

    vec.conservativeResize(mid_half_size);
}

void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::MatrixXd& matrix){
    const int size = pair_index.size();
    const int mid_half_size = size - 2*(size/4); // It is at least size/2
    std::vector<bool> index_bool(size, false);

    // Sort matrix columns according to pair_index
    for(int index=0; index<size;){
        int pos = pair_index[index].second;
        matrix.col(index).swap( matrix.col(pos) );
        index_bool[index] = true;

        int next = pair_index[pos].second;

        while(next != index){
            matrix.col(pos).swap(matrix.col(next));
            index_bool[pos] = true;

            pos = next;
            next = pair_index[pos].second;
        }

        index_bool[pos] = true;
        while(index_bool[index]) ++index;
    }

    // Move middle half to first part
    for(int index=0; index<mid_half_size;index++){
        matrix.col(index).swap(matrix.col(index+size/4));
    }

    matrix.conservativeResize(Eigen::NoChange , mid_half_size);
}

void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::MatrixXcd& matrix){
    const int size = pair_index.size();
    const int mid_half_size = size - 2*(size/4); // It is at least size/2
    std::vector<bool> index_bool(size, false);

    // Sort matrix columns according to pair_index
    for(int index=0; index<size;){
        int pos = pair_index[index].second;
        matrix.col(index).swap( matrix.col(pos) );
        index_bool[index] = true;

        int next = pair_index[pos].second;

        while(next != index){
            matrix.col(pos).swap(matrix.col(next));
            index_bool[pos] = true;

            pos = next;
            next = pair_index[pos].second;
        }

        index_bool[pos] = true;
        while(index_bool[index]) ++index;
    }

    // Move middle half to first part
    for(int index=0; index<mid_half_size;index++){
        matrix.col(index).swap(matrix.col(index+size/4));
    }
    matrix.resize(Eigen::NoChange , mid_half_size);
}
