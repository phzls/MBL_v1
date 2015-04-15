//
// Created by Liangsheng Zhang on 4/15/15.
//

#include <vector>
#include <Eigen/Dense>
#include "basis_transition.h"
#include "methods.h"

using namespace std;
using namespace Eigen;

void evec_to_basic(const EvolOP* evol, const vector<MatrixXd>& eigen, vector<vector<double> >& evec){
    if (evol -> Eigen_Basis_Type() == "Basic"){
        int index = 0;
        for (int i=0; i<eigen.size(); i++){
            if (eigen[i].rows() != evec[0].size()){
                cout << "The length of eigenvector in sector " << i << " is incompatible." << endl;
                cout << "Length from eigenvector: " << eigen[i].rows() << endl;
                cout << "Length from vector: " << evec[0].size() << endl;
                abort();
            }

            for (int j=0; j < eigen[i].cols(); j++){
                for (int k=0; k < eigen[i].rows(); k++){
                    evec[index][k] = eigen[i](k,j);
                }
                index ++;
            }
        }
    }

    else if (evol -> Eigen_Basis_Type() == "Parity"){
        // Assume only 2 sectors, 0 is even and 1 is odd
        // In the evec, even eigenvectors come first
        if (eigen.size() != 2){
            cout << "The position of even and odd sector in parity basis is not clear." << endl;
            abort();
        }

        TransitionMatrix transition;
        evol -> Transition_Compute(transition, "Basic_Parity", eigen);

        if (evec[0].size() != transition.Matrix("Basic_Even").rows()){
            cout << "The dimension of eigenvector acceptor is not correct." << endl;
            cout << "Assumed vector length: " << evec[0].size() << endl;
            cout << "Would obtained vector length: " << transition.Matrix("Basic_Even").rows()
            << endl;
            abort();
        }

        MatrixXd single_evec;

        int index = 0;
        for (int i=0; i< eigen[0].cols(); i++){
            // There may be a type problem, as transition matrix is always complex
            single_evec = transition.Matrix("Basic_Even") * eigen[0].col(i);

            for (int j=0; j< single_evec.size();j++)
                evec[index][j] = single_evec[j];

            index ++;
        }

        for (int i=0; i< eigen[1].cols(); i++){
            single_evec = transition.Matrix("Basic_Odd") * eigen[1].col(i);
            for (int j=0; j< single_evec.size();j++)
                evec[index][j] = single_evec[j];

            index ++;
        }
    }

    else{
        cout << "Eigen type " << evol -> Eigen_Basis_Type() << " unknown." << endl;
        abort();
    }
}

void evec_to_basic(const EvolOP* evol, const vector<MatrixXd>& eigen, vector<vector<complex<double> > >& evec){
    if (evol -> Eigen_Basis_Type() == "Basic"){
        int index = 0;
        for (int i=0; i<eigen.size(); i++){
            if (eigen[i].rows() != evec[0].size()){
                cout << "The length of eigenvector in sector " << i << " is incompatible." << endl;
                cout << "Length from eigenvector: " << eigen[i].rows() << endl;
                cout << "Length from vector: " << evec[0].size() << endl;
                abort();
            }

            for (int j=0; j < eigen[i].cols(); j++){
                for (int k=0; k < eigen[i].rows(); k++){
                    evec[index][k] = complex<double>(eigen[i](k,j),0);
                }
                index ++;
            }
        }
    }

    else if (evol -> Eigen_Basis_Type() == "Parity"){
        // Assume only 2 sectors, 0 is even and 1 is odd
        // In the evec, even eigenvectors come first
        if (eigen.size() != 2){
            cout << "The position of even and odd sector in parity basis is not clear." << endl;
            abort();
        }

        TransitionMatrix transition;
        evol -> Transition_Compute(transition, "Basic_Parity", eigen);

        if (evec[0].size() != transition.Matrix("Basic_Even").rows()){
            cout << "The dimension of eigenvector acceptor is not correct." << endl;
            cout << "Assumed vector length: " << evec[0].size() << endl;
            cout << "Would obtained vector length: " << transition.Matrix("Basic_Even").rows()
            << endl;
            abort();
        }

        MatrixXcd single_evec;

        int index = 0;
        for (int i=0; i< eigen[0].cols(); i++){
            // There may be a type problem, as transition matrix is always complex
            single_evec = transition.Matrix("Basic_Even") * eigen[0].col(i);

            for (int j=0; j< single_evec.size();j++)
                evec[index][j] = single_evec[j];

            index ++;
        }

        for (int i=0; i< eigen[1].cols(); i++){
            single_evec = transition.Matrix("Basic_Odd") * eigen[1].col(i);
            for (int j=0; j< single_evec.size();j++)
                evec[index][j] = single_evec[j];

            index ++;
        }
    }

    else{
        cout << "Eigen type " << evol -> Eigen_Basis_Type() << " unknown." << endl;
        abort();
    }
}

void evec_to_basic(const EvolOP* evol, const vector<MatrixXcd>& eigen, vector<vector<complex<double> > >& evec){
    if (evol -> Eigen_Basis_Type() == "Basic"){
        int index = 0;
        for (int i=0; i<eigen.size(); i++){
            if (eigen[i].rows() != evec[0].size()){
                cout << "The length of eigenvector in sector " << i << " is incompatible." << endl;
                cout << "Length from eigenvector: " << eigen[i].rows() << endl;
                cout << "Length from vector: " << evec[0].size() << endl;
                abort();
            }

            for (int j=0; j < eigen[i].cols(); j++){
                for (int k=0; k < eigen[i].rows(); k++){
                    evec[index][k] = eigen[i](k,j);
                }
                index ++;
            }
        }
    }

    else if (evol -> Eigen_Basis_Type() == "Parity"){
        // Assume only 2 sectors, 0 is even and 1 is odd
        // In the evec, even eigenvectors come first
        if (eigen.size() != 2){
            cout << "The position of even and odd sector in parity basis is not clear." << endl;
            abort();
        }

        TransitionMatrix transition;
        evol -> Transition_Compute(transition, "Basic_Parity", eigen);

        if (evec[0].size() != transition.Matrix("Basic_Even").rows()){
            cout << "The dimension of eigenvector acceptor is not correct." << endl;
            cout << "Assumed vector length: " << evec[0].size() << endl;
            cout << "Would obtained vector length: " << transition.Matrix("Basic_Even").rows()
            << endl;
            abort();
        }

        MatrixXcd single_evec;

        int index = 0;
        for (int i=0; i< eigen[0].cols(); i++){
            // There may be a type problem, as transition matrix is always complex
            single_evec = transition.Matrix("Basic_Even") * eigen[0].col(i);

            for (int j=0; j< single_evec.size();j++)
                evec[index][j] = single_evec[j];

            index ++;
        }

        for (int i=0; i< eigen[1].cols(); i++){
            single_evec = transition.Matrix("Basic_Odd") * eigen[1].col(i);
            for (int j=0; j< single_evec.size();j++)
                evec[index][j] = single_evec[j];

            index ++;
        }
    }

    else{
        cout << "Eigen type " << evol -> Eigen_Basis_Type() << " unknown." << endl;
        abort();
    }
}

