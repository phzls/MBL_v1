//
// Created by Liangsheng Zhang on 4/14/15.
//

#ifndef MBL_V1_SCREEN_OUTPUT_H
#define MBL_V1_SCREEN_OUTPUT_H

/*
 * This file includes some functions which print out contents to screen.
 */

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <cmath>

using namespace std;
using namespace Eigen;


/*
 * Write out a complex number. If the imaginary part has magnitude smaller than delta, then
 * it is assumed to 0 and not written out. If the imaginary part has magnitude smaller than
 * delta, then it is assumed to 0 and not written out. If both real and imaginary parts are
 * not written out, then a 0 is outputted.
 */
template <class T>
void complex_write(const complex<T>& val){
    const double delta = 1.0e-15;

    T re = real(val);
    bool re_out = true;

    T im = imag(val);
    bool im_out = true;

    if (abs(re) > delta){
        cout << real(val);
    }
    else re_out = false;

    if(abs(im)>delta){
        if (im>=0 && re_out) cout << "+";
        cout << im << "j";
    }
    else im_out = false;

    if(!re_out && !im_out) cout << 0;

    cout << "    ";
}

/*
 * Return a complex number in string. If the imaginary part has magnitude smaller than delta, then
 * it is assumed to 0 and not written out. If the imaginary part has magnitude smaller than
 * delta, then it is assumed to 0 and not written out. If both real and imaginary parts are
 * not written out, then a 0 is outputted.
 */
template <class T>
string complex_write_return(const complex<T>& val){
    stringstream c_num;
    const double delta = 1.0e-15;

    T re = real(val);
    bool re_out = true;

    T im = imag(val);
    bool im_out = true;

    if (abs(re) > delta){
        c_num << real(val);
    }
    else re_out = false;

    if(abs(im)>delta){
        if (im>=0 && re_out) c_num << "+";
        c_num << im << "j";
    }
    else im_out = false;

    if(!re_out && !im_out) c_num << 0;

    return c_num.str();
}

/*
 * Use complex_write function to write out an Eigen dynamic complex matrix
 */
void complex_matrix_write(const MatrixXcd&);

/*
 * Write out an Eigen dynamic real matrix
 */
void real_matrix_write(const MatrixXd&);

/*
 * Write out an Eigen dynamic matrix
 */
void matrix_write(const MatrixXcd&);
void matrix_write(const MatrixXd&);


#endif //MBL_V1_SCREEN_OUTPUT_H
