//
// Created by Liangsheng Zhang on 9/18/15.
//

#include <iostream>
#include <vector>
#include "disorder_transition_func.h"
#include "screen_output.h"

using namespace std;

void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::VectorXd& vector){
    const int size = pair_index.size();
    const int mid_half_size = size - 2*(size/4); // It is at least size/2
    std::vector<bool> index_bool(mid_half_size, false);

    int step = 0;
    for(int index=0; index<mid_half_size;){
        int pos = pair_index[(index + size/4)%size].second;

        swap(vector(index), vector(pos));
        index_bool[index] = true;

        int next = pair_index[(pos+size/4)%size].second;
        int init_obj_pos = pos;

        while(next != index){
            if(pos<mid_half_size){
                swap(vector(pos), vector(next));
                init_obj_pos = next;
                index_bool[pos] = true;
            }
            else if(next < mid_half_size){
                swap(vector(next), vector(init_obj_pos));
                init_obj_pos = next;
                index_bool[next] = true;
            }

            pos = next;
            next = pair_index[(pos+size/4)%size].second;
        }

        while(index_bool[index]) ++index;
    }

    vector.conservativeResize(mid_half_size);
}

void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::MatrixXd& matrix){
    const int size = pair_index.size();
    const int mid_half_size = size - 2*(size/4); // It is at least size/2
    std::vector<bool> index_bool(mid_half_size, false);

    for(int index=0; index<mid_half_size;){
        int pos = pair_index[(index + size/4)%size].second;
        matrix.col(index).swap( matrix.col(pos) );
        index_bool[index] = true;

        int next = pair_index[(pos+size/4)%size].second;
        int init_obj_pos = pos;

        while(next != index){
            if(pos<mid_half_size){
                matrix.col(pos).swap(matrix.col(next));
                init_obj_pos = next;
                index_bool[pos] = true;
            }
            else if(next < mid_half_size){
                matrix.col(next).swap( matrix.col(init_obj_pos) );
                init_obj_pos = next;
                index_bool[next] = true;
            }

            pos = next;
            next = pair_index[(pos+size/4)%size].second;
        }

        while(index_bool[index]) ++index;
    }

    matrix.conservativeResize(Eigen::NoChange , mid_half_size);
}

void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::MatrixXcd& matrix){
    const int size = pair_index.size();
    const int mid_half_size = size - 2*(size/4); // It is at least size/2
    std::vector<bool> index_bool(mid_half_size, false);

    for(int index=0; index<mid_half_size;){
        int pos = pair_index[(index + size/4)%size].second;
        matrix.col(index).swap( matrix.col(pos) );
        index_bool[index] = true;

        int next = pair_index[(pos+size/4)%size].second;
        int init_obj_pos = pos;

        while(next != index){
            if(pos<mid_half_size){
                matrix.col(pos).swap(matrix.col(next));
                init_obj_pos = next;
                index_bool[pos] = true;
            }
            else if(next < mid_half_size){
                matrix.col(next).swap( matrix.col(init_obj_pos) );
                init_obj_pos = next;
                index_bool[next] = true;
            }

            pos = next;
            next = pair_index[(pos+size/4)%size].second;
        }

        while(index_bool[index]) ++index;
    }

    matrix.resize(Eigen::NoChange , mid_half_size);
}
