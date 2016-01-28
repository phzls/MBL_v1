//
// Created by Liangsheng Zhang on 1/27/16.
//

/*
 * This function computes entropies of all eigenstates from a spin chain model
 * where the local dimension is 2. The entropy concerns about the left part
 * of the spin which has left_size. The eigenvectors can be either real or
 * imaginary and they are converted to basic binary basis representation
 * for computation of reduced density matrix.
 */

#include <iostream>
#include "methods.h"

using namespace std;
using namespace Eigen;

void model_entropy_left_2(const EvolOP* model, int left_size, vector<double>& ent){
    const int dim = model->Get_Dim();
    const int size = model->Get_Size();

    if( (1<<size) != dim){
        cout << "model " << model->Repr() << " may not have local dimension 2" << endl;
        cout << "size: " << size << " dim: " << dim << " Expected dim: " << (1<<dim) << endl;
        abort();
    }

    ent.resize(dim);
    for(int i=0; i<ent.size(); i++) ent[i] = 0;

    if(model->Evec_Type() == "Real"){
        // Vector for eigenvectors in basic binary basis
        vector<vector<double > > evec_basic(dim);
        for (int i=0; i<evec_basic.size();i++) evec_basic[i].resize(dim);

        vector<MatrixXd>* evec_real = new vector<MatrixXd>;
        model->Evec(*evec_real);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, *evec_real, evec_basic);

        // Release memory
        delete evec_real;

        VectorXd evec(dim); // Vector for one eigenstate in binary basis
        for(int i=0; i<evec_basic.size(); i++){
            for(int j=0; j<evec_basic[i].size();j++) evec[j] = evec_basic[i][j];
            ent[i] = state_entropy_left_2(evec, size, left_size);
        }
    }
    else if(model->Evec_Type() == "Complex"){
        // Vector for eigenvectors in basic binary basis
        vector<vector< complex<double> > > evec_basic(dim);
        for (int i=0; i<evec_basic.size();i++) evec_basic[i].resize(dim);

        vector<MatrixXcd>* evec_complex = new vector<MatrixXcd>;
        model->Evec(*evec_complex);

        // Convert eigenvectors in basic basis
        evec_to_basic(model, *evec_complex, evec_basic);

        // Release memory
        delete evec_complex;

        VectorXcd evec(dim); // Vector for one eigenstate in binary basis
        for(int i=0; i<evec_basic.size(); i++){
            for(int j=0; j<evec_basic[i].size();j++) evec[j] = evec_basic[i][j];
            ent[i] = state_entropy_left_2(evec, size, left_size);
        }
    }
    else{
        cout << "Eigenvector type " << model->Evec_Type() << " is not recognized" << endl;
        abort();
    }
}

