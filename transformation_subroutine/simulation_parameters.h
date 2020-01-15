//
// Created by erolsson on 25/09/2019.
//

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include "Eigen/Dense"

class TransformationMaterialParameters {
public:
    TransformationMaterialParameters(const double* data, double hv):
            data_(data), back_stresses_(static_cast<int>(data[10])), hv_(hv) {}

    // Elastic parameters
    [[nodiscard]] const double& E() const { return data_[0]; }         // Young's modulus
    [[nodiscard]] const double& v() const { return data_[1]; }         // Poisson's ratio

    // Initial yield stress of Martensite and Austenite
    [[nodiscard]] double sy0M() const { return data_[2] + data_[3]*hv_; } ;     // Initial yield stress
    [[nodiscard]] double sy0A() const { return data_[4] + data_[5]*hv_; } ;     // Initial yield stress

    // Parameters for isostropic hardnening
    [[nodiscard]] double Q() const { return data_[6] + data_[7]*hv_; } ;           // Maximum size of yield surface
    [[nodiscard]] double b() const { return data_[8] + data_[9]*hv_; } ;           // Growth rate of yield surface

    // Parameters for kinematic hardnening
    [[nodiscard]] unsigned back_stresses() const { return back_stresses_; }
    [[nodiscard]] double Cm(unsigned n) const { return data_[11 + 4*n] + data_[11 + 4*n + 1]*hv_; }
    [[nodiscard]] double gamma(unsigned n) const { return data_[11 + 4*n + 2] + data_[11 + 4*n + 3]*hv_; }

    // Parameter for the SDE effect
    [[nodiscard]] const double& a() const { return data_[11 + 4*back_stresses_]; } ;

    // Parameters for Greenwood Johnson
    [[nodiscard]] const double& R1() const { return data_[12 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& R2() const { return data_[13 + 4*back_stresses_]; } ;

    // Volume expansion of martensite compared to austenite
    [[nodiscard]] const double& dV() const { return data_[14 + 4*back_stresses_]; } ;

    // Martensite start temperature
    [[nodiscard]] const double& Ms() const { return data_[15 + 4*back_stresses_]; } ;

    // Stabilisation temperature due to tempering
    [[nodiscard]] const double& Mss() const { return data_[16 + 4*back_stresses_]; } ;

    // Koistinen-Marburger parameter
    [[nodiscard]] const double& k() const { return data_[17 + 4*back_stresses_]; } ;

    // Material parameters for stress induced phase transformations
    [[nodiscard]] const double& a1() const { return data_[18 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& a2() const { return data_[19 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& a3() const { return data_[20 + 4*back_stresses_]; } ;

    // Material parameters for the strain induced phase transformations
    [[nodiscard]] const double& beta() const { return data_[21 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& alpha() const { return data_[22 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& n() const { return data_[23 + 4*back_stresses_]; } ;

    [[nodiscard]] const double& g0() const { return data_[24 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& g1() const { return data_[25 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& g2() const { return data_[26 + 4*back_stresses_]; } ;

    [[nodiscard]] const double& g_mean() const { return data_[27 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& g_std() const { return data_[28 + 4*back_stresses_]; } ;

    [[nodiscard]] const double& M_sigma() const { return data_[29 + 4*back_stresses_]; } ;
    [[nodiscard]] const double& M_d() const { return data_[30 + 4*back_stresses_]; } ;

    [[nodiscard]] bool plastic() const { return sy0M() > 0 && sy0A() > 0; }
    [[nodiscard]] bool kinematic_hardening() const { return back_stresses_ > 0; }

private:
    const double* data_;
    unsigned back_stresses_ = 0;
    double hv_;
};


#endif // SIMULATION_PARAMETERS_H
