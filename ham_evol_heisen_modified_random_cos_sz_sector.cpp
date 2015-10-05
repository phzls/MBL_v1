//
// Created by Liangsheng Zhang on 10/5/15.
//

/**
 ** Implementation of HamEvolHeisenModifiedRandomCosSzSector class
 **/

#include <complex>
#include <iostream>
#include <algorithm>
#include <map>
#include "constants.h"
#include "ham_evol_model.h"
#include "screen_output.h"
#include "randomc.h"
#include "combinatorics.h"

using namespace std;

extern CRandomMersenne RanGen_mersenne; // points in [0,1)

/*
 * Construct the representation string and abstract type of the class.
 */
void HamEvolHeisenModifiedRandomCosSzSector::Repr_Init_(){
    repr_ << "Heisen_Modified_Random_Cos_Sz_Sector_Hamiltonian_L=" << size_ << ",h=" << h_ << ",Sz_" << total_spin_z_;
    type_ = "Heisen_Modified_Random_Cos_Sz_Sector_Hamiltonian";
}

/*
 * Initialize random numbers
 */
void HamEvolHeisenModifiedRandomCosSzSector::Evol_Para_Init() {
    random_h_.resize(size_);
    vector<double> theta(size_);

    // Random Fields
    for (int i=0; i<size_;i++){
        // Random phase
        double u = RanGen_mersenne.Random();
        theta[i] = 2 * Pi * u;
        random_h_[i] = h_ * cos(theta[i]);
    }

    if (debug_){
        cout << "Random phase: " << endl;
        for (int i=0; i<size_; i++) cout << theta[i] << endl;
        cout << "Longitude fields:" << endl;
        for (int i=0; i<size_; i++) cout << random_h_[i] << endl;
        cout << endl;
    }
}

/*
 * Compute the dimension of the restricted Hilbert space for that sector
 */
void HamEvolHeisenModifiedRandomCosSzSector::Dim_Cal_() {
    int temp = total_spin_z_ + size_;
    if( temp % 2 != 0){
        cout << "Cannot achieve total Z spin " << total_spin_z_ << " for system size " << size_ << endl;
        abort();
    }

    int up_spin = temp/2; // Number of up spins
    dim_ = binomial_coef(size_, up_spin);
}

/*
 * Construct availabe spin configurations: it gives basic basis numbers whose total
 * Z spin agrees with the specified number
 * The config_ vector will be kept even when the Hamiltonian is destroyed to used
 * for possible future basis transformation
 */

void HamEvolHeisenModifiedRandomCosSzSector::Spin_Config_() {
    if(!config_constructed_) { // Otherwise no need to construct
        int temp = total_spin_z_ + size_;
        if( temp % 2 != 0){
            cout << "Cannot achieve total Z spin " << total_spin_z_ << " for system size " << size_ << endl;
            abort();
        }

        int up_spin = temp/2; // Number of up spins
        comb_num(0, size_-1, up_spin, config_);

        if(config_.size() != dim_){
            cout << "Incompatible dimension for " << repr_ << endl;
            cout << "Number of configurations: " << config_.size() << endl;
            cout << "Dimension of restricted Hilbert space: " << dim_ << endl;
            abort();
        }

        config_constructed_ = true;

        // Sort in ascending order
        sort(config_.begin(), config_.end());

        if(debug_){
            cout << "Spin states in basic basis with total z spin " << total_spin_z_ << ":" << endl;
            for(int i=0; i<config_.size(); i++)
                cout << config_[i] << endl;
        }
    }
}

/*
 * Construct the Hamiltonian only for the sector of given total Z spin
 */
void HamEvolHeisenModifiedRandomCosSzSector::Evol_Construct() {

    if (!constructed_){
        ham_op_ = MatrixXd::Zero(dim_,dim_);
        constructed_ = true;
    }

    else{
        cout << "Evolution operator has been constructed." << endl;
        abort();
    }
    if(size_ < 1) {
        cout << "Cannot construct Hamiltonian with " << size_ << " number of spins for " << type_ << endl;
        abort();
    }

    Spin_Config_();

    // Change config_ to a map for easy lookup of basic states to restricted states
    map<int, int> config_map;
    for(int i=0; i<config_.size(); i++) {
        int state = config_[i];
        config_map[state] = i;
    }

    for(int i=0; i<config_.size(); i++){
        int state = config_[i];
        int mod_state = state; // The state that would be modified later
        int spin = 0;
        int prev_spin = 0;
        int flip_state = -1; // The state after some spins are flipped in full basic basis
        int restricted_flip_state = -1; // Corresponding state in restricted basic basis
        // for flip_state

        int pos = 1; // Representing the position of spin
        int prev_pos = 1; // Representing the position of prev_spin

        double z_val = 0; // diagonal term

        // The leftmost spin
        spin = 2*(mod_state & 1) - 1;
        mod_state = mod_state >> 1;
        z_val += random_h_[0] * spin;

        for(int j=1; j<size_; j++){
            prev_spin = spin;
            prev_pos = pos;

            spin = 2*(mod_state & 1) - 1;
            mod_state = mod_state >> 1;
            pos = pos << 1;

            z_val += random_h_[j] * spin;
            z_val += spin * prev_spin;

            if(prev_spin + spin == 0){
                // We should stay in the same sector of given total Z spin
                // so the flipped two spins must point to opposite directions
                flip_state = state ^ prev_pos; // flip prev_spin
                flip_state = flip_state ^ pos; // flip spin

                restricted_flip_state = config_map[flip_state];
                ham_op_(restricted_flip_state, i) += 2*(1-h_);
            }
            ham_op_(i, i) = z_val;
        }
    }

    // Check if the matrix is Hermitian
    for (int i=0; i<ham_op_.cols(); i++){
        for (int j=i+1; j< ham_op_.rows();j++){
            if (abs(ham_op_(j,i) - ham_op_(i,j)) > 1.0e-10){
                cout << "Hamiltonian is not Hermitian at row " << j <<" and col " << i << endl;
                cout << "(i,j): " << ham_op_(i,j) << endl;
                cout << "(j,i): " << ham_op_(j,i) << endl;
                abort();
            }
        }
    }

    if (debug_){
        cout << "Hamiltonian matrix:" << endl;
        matrix_write(ham_op_);
        cout << endl;
    }
}

