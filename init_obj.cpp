//
// Created by Liangsheng Zhang on 6/1/15.
//

#include <iostream>
#include <map>
#include <sstream>
#include "init_obj.h"

using namespace std;

InitInfo::InitInfo(const InitInfo& init_info){
    Copy(init_info);
}

void InitInfo::Copy(const InitInfo & init_info) {
    size = init_info.size;
    dim = init_info.dim;
    norm_delta = init_info.norm_delta;
    debug = init_info.debug;
    init_func_name = init_info.init_func_name;
}

void InitObj::Init_Func(const TransitionMatrix& transition, VectorXcd& init_state) const {

    if (init_info.init_func_name == ""){
        cout << "Function name for initial state construction has not been initialized." << endl;
        abort();
    }

    map<string, init_func>::const_iterator it = init_func_map_.find(init_info.init_func_name);

    if (it == init_func_map_.end()){
        cout << "The requested init_func does not exist." << endl;
        Print();
        cout << "Requested function: " << init_info.init_func_name << endl;
        abort();
    }
    else it -> second(init_info, transition, init_state);
}

void InitObj::Init_Func_C(const TransitionMatrix& transition, MatrixXcd& density) const {

    if (init_info.init_func_name == ""){
        cout << "Function name for initial density matrix construction has not been initialized." << endl;
        abort();
    }

    map<string, init_func_C>::const_iterator it = init_func_C_map_.find(init_info.init_func_name);

    if (it == init_func_C_map_.end()){
        cout << "The requested init_func does not exist." << endl;
        Print_C();
        cout << "Requested function: " << init_info.init_func_name << endl;
        abort();
    }
    else it -> second(init_info, transition, density);
}

void InitObj::Print() const {
    map<string, init_func>::const_iterator it;

    cout << "Initial State Construction Functions: "<< endl;
    for (it = init_func_map_.begin(); it != init_func_map_.end(); it++){
        cout << it -> first << endl;
    }
}

void InitObj::Print_C() const {
    map<string, init_func_C>::const_iterator it;

    cout << "Initial Density Construction Functions: "<< endl;
    for (it = init_func_C_map_.begin(); it != init_func_C_map_.end(); it++){
        cout << it -> first << endl;
    }
}

void InitObj::map_init_(){
    map<string, init_func>::const_iterator it;
    map<string, init_func_C>::const_iterator it_C;

    string name1 = "Random_Product";
    init_func func1 = random_product;
    init_func_C func_C1 = random_product;

    it = init_func_map_.find(name1);
    if (it != init_func_map_.end()){
        cout << "init_func " << name1 << " already exists." << endl;
        abort();
    }
    init_func_map_[name1] = func1;

    it_C = init_func_C_map_.find(name1);
    if (it_C != init_func_C_map_.end()){
        cout << "init_func " << name1 << " already exists." << endl;
        abort();
    }
    init_func_C_map_[name1] = func_C1;

    string name2 = "Product_Random";
    init_func func2 = product_random;

    it = init_func_map_.find(name2);
    if (it != init_func_map_.end()){
        cout << "init_func " << name2 << " already exists." << endl;
        abort();
    }
    init_func_map_[name2] = func2;

    string name3 = "Random_Pure";
    init_func_C func_C3 = random_pure;

    it_C = init_func_C_map_.find(name3);
    if (it_C != init_func_C_map_.end()){
        cout << "init_func " << name3 << " already exists." << endl;
        abort();
    }
    init_func_C_map_[name3] = func_C3;

    string name4 = "Left_Spin_Random";
    init_func_C func_C4 = left_spin_random;

    it_C = init_func_C_map_.find(name4);
    if (it_C != init_func_C_map_.end()){
        cout << "init_func " << name4 << " already exists." << endl;
        abort();
    }
    init_func_C_map_[name4] = func_C4;
}