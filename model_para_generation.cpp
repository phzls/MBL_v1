//
// Created by Liangsheng Zhang on 1/27/16.
//

/*
 * This function generates random parameters used for various models
 */

#include <iostream>
#include <cmath>
#include "constants.h"
#include "methods.h"
#include "randomc.h"


using namespace std;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

void model_para_generation(const AllPara& parameters, vector< vector<double> >& model_para, int size){
    const string model_name = parameters.generic.model;
    const bool debug = parameters.generic.debug;

    if(model_name == "XXZ_Uniform_Z_Random_Shift_Real_Flo"){
        model_para.resize(1);
        model_para[0].resize(size);

        for (int i=0; i<size; i++){
            double u = RanGen_mersenne.Random();
            model_para[0][i] = 2*sqrt(3)*u - sqrt(3);
        }
    }
    else if(model_name == "XXZ_Gaussian_Z_Random_Shift_Real_Flo"){
        model_para.resize(1);
        model_para[0].resize(size);

        for (int i=0; i<size;){
            // Using Box Muller transformation to general standard normal
            double u1 = RanGen_mersenne.Random();
            double u2 = RanGen_mersenne.Random();
            model_para[0][i] = sqrt(-2*log(u1))*cos(2*Pi*u2);
            i++;
            if(i<size) model_para[0][i] = sqrt(-2*log(u1))*sin(2*Pi*u2);
            i++;
        }
    }
    else{
        cout << "Paramter generation for model " << model_name << " is not implemented yet" << endl;
        abort();
    }

    if(debug){
        cout << "For model " << model_name << " random numbers:" << endl;
        for(int i=0; i<model_para.size(); i++){
            for(int j=0; j<model_para[i].size(); j++) cout << model_para[i][j] << endl;
            cout << endl;
        }
    }
}

