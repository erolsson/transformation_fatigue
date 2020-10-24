//
// Created by erolsson on 11/09/2019.
//

#ifndef UTILITIES_H
#define UTILITIES_H

#include "Eigen/Dense"

Eigen::IOFormat CleanFmt(8, 0, ", ", "\n", "[", "]");

using Eigen::Matrix;
// Unit forth order tensor in matrix form
const static Matrix<double, 6, 6> I((Matrix<double, 6, 6>() <<
        1.,  0.,  0.,  0.,  0.,  0.,
        0.,  1.,  0.,  0.,  0.,  0.,
        0.,  0.,  1.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.5, 0.,  0.,
        0.,  0.,  0.,  0.,  0.5, 0.,
        0.,  0.,  0.,  0.,  0.,  0.5).finished());

//Deviatoric fourth order tensor in matrix form
const static Matrix<double, 6, 6> J((Matrix<double, 6, 6>() <<
         2./3,  -1./3,  -1./3,  0.,  0.,  0.,
        -1./3.,  2./3,  -1./3,  0.,  0.,  0.,
        -1./3,  -1./3,   2./3., 0.,  0.,  0.,
         0.,     0.,     0.,    0.5, 0.,  0.,
         0.,     0.,     0.,    0.,  0.5, 0.,
         0.,     0.,     0.,    0.,  0.,  0.5).finished());

const static Matrix<double, 6, 6> E3((Matrix<double, 6, 6>() <<
        1.,  1.,  1.,  0.,  0.,  0.,
        1.,  1.,  1.,  0.,  0.,  0.,
        1.,  1.,  1.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.).finished());

// Kronecker delta in vector notation
const static Matrix<double, 6, 1> delta_ij((Matrix<double, 6, 1>() <<
       1.,  1.,  1.,  0.,  0.,  0.).finished());

template<typename T>
Matrix<T, 6, 1> deviator(const Matrix<T, 6, 1>& tensor) {
    T hydrostatic = (tensor[0] + tensor[1] + tensor [2])/3;
    Matrix<T, 6, 1> dev = tensor;
    dev[0] -= hydrostatic;
    dev[1] -= hydrostatic;
    dev[2] -= hydrostatic;
    return dev;
}

// Double contraction for the operation a_ij b_ij
template<typename T>
T double_contract(const Matrix<T, 6, 1>& a, const Matrix<T, 6, 1>& b) {
    T val = 0;
    for (unsigned i = 0; i!= 6; ++i) {
        val += a[i]*b[i];
    }
    val += a[3]*b[3] + a[4]*b[4] + a[5]*b[5];
    return val;
}

// Double contraction for the operation A_ijkl b_kl
template<typename T>
Matrix<T, 6, 1> double_contract(const Matrix<T, 6, 6>& A,  const Matrix<T, 6, 1>& b) {
    Matrix<T, 6, 1> result = A*b;
    for (unsigned i = 0; i!= 6; ++i) {
        for (unsigned j = 3; j != 6; ++j) {
            result[i] += A(i, j)*b[j];
        }
    }
    return result;
}

// Double contraction for the operation A_ijkl B_klmn
template<typename T>
Matrix<T, 6, 6> double_contract(const Matrix<T, 6, 6>& A,
                                                 const Matrix<T, 6, 6>& B) {
    Matrix<T, 6, 6> res = A*B;
    for (unsigned i = 0; i != 6; ++i) {
        for (unsigned j = 0; j != 6; ++j) {
            for (unsigned k = 3; k != 6; ++k) {
                res(i, j) += A(i, k)*B(k, j);
            }
        }
    }
    return res;
}

// contraction for the operation a_ij b_jk
template<typename T>
Matrix<T, 6, 1> contract(const Matrix<T, 6, 1> a, const Matrix<T, 6, 1>& b) {
    Matrix<T, 6, 1> result = Matrix<T, 6, 1>::Zero();
    result[0] = a[0]*b[0] + a[3]*b[3] + a[4]*b[4];
    result[1] = a[1]*b[1] + a[3]*b[3] + a[5]*b[5];
    result[2] = a[2]*b[2] + a[4]*b[4] + a[5]*b[5];
    result[3] = a[0]*b[3] + a[3]*b[1] + a[4]*b[5];
    result[4] = a[0]*b[4] + a[3]*b[5] + a[4]*b[2];
    result[5] = a[3]*b[4] + a[1]*b[5] + a[5]*b[2];
    return result;
}

template<typename T>
T von_Mises(const Matrix<T, 6, 1>& s) {
    double vm2 = s[0]*s[0] + s[1]*s[1] + s[2]*s[2] - s[0]*s[1] - s[0]*s[2] - s[1]*s[2] +
                 3*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
    // vm2 could be a very small negative number due to round off errors
    if (vm2 > 0) {
        return sqrt(vm2);
    }
    return 0.;
}

template<typename T>
T vector_det(const Matrix<T, 6, 1>& s) {
    return s[0]*(s[1]*s[2] - s[5]*s[5]) - s[3]*(s[3]*s[2] - s[5]*s[4]) + s[4]*(s[3]*s[5] - s[1]*s[4]);
}

#endif //UTILITIES_H
