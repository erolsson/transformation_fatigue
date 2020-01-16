//
// Created by erolsson on 11/09/2019.
//

#include "transformation_umat.h"
#include "transformation_subroutine_header.h"

#include <cfenv>
#include "fenv.h"
#include <cmath>
#include <iostream>

#include "Eigen/Dense"

#include "simulation_parameters.h"
#include "stress_functions.h"

#include <mutex>

#pragma STDC FENV_ACCESS ON
const double pi = 3.14159265359;
std::mutex print_mutex;

/*
void print_at_time(const std::string msg, double time, unsigned noel, unsigned npt) {
    if (time > -1) {
        std::cout << "Elem: " << noel << " gp: " << npt << " " << msg << std::endl;
    }
}
*/

template<typename T>
void print_for_position(const std::string msg, const T& val, unsigned noel, unsigned npt) {
    //if (noel == 57059 && npt == 5) {
    std::cout << "noel=" << noel << "  gp=" <<  npt << "  " << msg << val << std::endl;
    //}
}

class State {
public:
    using Vector6 = Eigen::Matrix<double, 6, 1>;
    explicit State(double* data, unsigned back_stresses, double* drot) :
        data_(data), back_stresses_(back_stresses) {
        for(unsigned i = 0; i != back_stresses_; ++i) {
            double alpha[6] = {0, 0, 0, 0, 0, 0};
            rotsig_(&data_[15 + i*6], drot, alpha, 1, 3, 3);
            for(unsigned j = 0; j != 6; ++j) {
                data_[15 + i*6 + j] = alpha[j];
            }
        }
    }
    double& ep_eff() { return data_ [0]; }
    double& austenite() { return  data_[1]; }
    double& carbon() { return  data_[2]; }

    double& ferrite() { return  data_[3]; }
    [[nodiscard]] const double& ferrite() const  { return  data_[3]; }

    double& hardness_hv() { return  data_[4]; }

    double& lbainite() { return  data_[5]; }
    [[nodiscard]] const double& lbainite() const  { return  data_[5]; }

    double& pearlite() { return  data_[6]; }
    [[nodiscard]] const double& pearlite() const  { return  data_[6]; }

    double& q_martensite() { return  data_[7]; }
    [[nodiscard]] const double& q_martensite() const { return  data_[7]; }

    double& t_martensite() { return  data_[8]; }
    [[nodiscard]] const double& t_martensite() const { return  data_[8]; }

    double& ubainite() { return  data_[9]; }
    [[nodiscard]] const double& ubainite() const { return  data_[9]; }

    [[nodiscard]] double fM()  const { return q_martensite() + t_martensite(); }
    [[nodiscard]] double other_phases() const { return ferrite() + lbainite() + pearlite() + ubainite(); }
    double& fM_stress() { return data_[10]; }
    double& fM_strain() { return data_[11]; }
    double& fsb() { return data_[12]; }
    double& fsb0() { return data_[13]; }
    double& R() { return data_ [14]; }

    Eigen::Map<Vector6> back_stress_vector(unsigned n) {
        return Eigen::Map<Vector6>(&data_[15 + n*6]);
    }

    Eigen::Map<Vector6> total_back_stress() {
        return Eigen::Map<Vector6>(&data_[15 + back_stresses_*6]);
    }

private:
    double* data_;
    unsigned back_stresses_;
};

double yield_function(const Eigen::Matrix<double, 6, 1>& sigma, const Eigen::Matrix<double, 6, 1>& alpha,
        double sigma_y, const TransformationMaterialParameters& params) {
    Eigen::Matrix<double, 6, 1> stilde = deviator(sigma) - alpha;
    double I1 = sigma[0] + sigma[1] + sigma[2];
    return sqrt(3*double_contract(stilde, stilde)/2) + params.a()*I1 - sigma_y;
}

double stress_temperature_transformation(const Eigen::Matrix<double, 6, 1>& stress,
                                         const TransformationMaterialParameters& params, double T) {
    Eigen::Matrix<double, 6, 1> s_dev = deviator(stress);
    double m_stress = params.a1()*(stress[0] + stress[1] + stress[2]);   // Contribution from hydrostatic stress
    m_stress += params.a2()*von_Mises(stress);
    m_stress += params.a3()*vector_det(s_dev);
    return params.k()*(params.Ms() + m_stress + params.Mss() - T);
}

double stress_transformation_function(const Eigen::Matrix<double, 6, 1>& stress, double T,
                                      const TransformationMaterialParameters& params, const State& state, double fM) {
    return (1 - exp(-stress_temperature_transformation(stress, params, T)))*(1 - state.other_phases()) - fM;
}

double normal_pdf(double x) {
    return exp(-0.5*x*x)/sqrt(2*pi);
}

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
        double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran, double* time,
        double& dtime, double& temp, double& dtemp, double *predef, double *dpred, char *cmname, const int& ndi,
        const int& nshr, const int& ntens, const int& nstatv, const double* props, const int& nprops, double *coords,
        double* drot, double& pnewdt, double& celent, double* dfgrd0, double* dfgrd1, const int& noel, const int& npt,
        const int& layer, const int& kspt, const int& kstep, const int& kinc, short cmname_len) {

    feenableexcept(FE_INVALID | FE_OVERFLOW);
    // std::lock_guard<std::mutex> lock(print_mutex);
    // print_at_time("starting", time[1], noel, npt);
    using Matrix6x6 = Eigen::Matrix<double, 6, 6>;
    using Vector6 = Eigen::Matrix<double, 6, 1>;

    const TransformationMaterialParameters params(props, statev[4]);
    // print_at_time("setting params", time[1], noel, npt);
    State state(statev, params.back_stresses(), drot);
    // print_at_time("setting state", time[1], noel, npt);
    Eigen::Map<Vector6> stress_vec(stress);

    const Eigen::Map<Vector6> de(dstran);
    // print_for_position("strain_inc: ", de.transpose().format(CleanFmt), noel, npt);
    Eigen::Map<Matrix6x6> D_alg(ddsdde);

    // Elastic parameters
    double G = params.E()/2/(1 + params.v());
    double K = params.E()/3/(1 - 2*params.v());
    Matrix6x6 Del = 2*G*J + K*E3;

    // Yield stress at start of increment
    double sy = params.sy0M()*state.fM() + params.sy0A()*(1 - state.fM()) + state.R();
    Vector6 sigma_t = stress_vec + Del*de;  // Trial stress
    Vector6 sij_t = deviator(sigma_t);

    Vector6 stilde2 = sij_t;
    if (params.kinematic_hardening()) {
        stilde2 -= state.total_back_stress();
    }
    bool plastic = params.plastic() && yield_function(sigma_t, state.total_back_stress(), sy, params) > 0;
    // print_at_time("yield function evaluated", time[1], noel, npt);
    bool stress_transformations = stress_transformation_function(sigma_t, temp, params, state, state.fM()) >= 0;
    bool strain_transformations = params.beta() > 0 && plastic;
    bool elastic = !plastic && !stress_transformations;

    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        // print_at_time("Elastic increment", time[1], noel, npt);
        D_alg = Del;
        stress_vec = sigma_t;
    }
    else {  // Inelastic deformations
        // print_at_time("Non elastic increment", time[1], noel, npt);
        // Increment in plastic strain and martensitic phase fraction
        Vector6 sigma_2 = sigma_t;
        Vector6 s = deviator(sigma_2);
        double DL = 0;
        double DfM_stress = 0;
        double DfM_strain = 0;

        double DfM = 0;

        double f = 0;
        double h_stress = 0;
        double h_strain = 0;

        double dfdDL = 0;
        double dfdDfM = 0;
        double dh_stressDL = 0;
        double dh_stressDfM = 0;

        double dh_straindDL = 0;
        double dh_strainDfM = 0;

        double fsb2 = state.fsb();

        double As = 0;
        double Bs = 0;
        double DSigma = 0;
        double Sigma = 0;
        double P = 0;
        double I1_2 = 0;
        double pdf = 0;
        double dfsb2dDL = 0;
        double norm_drivning_force = 0;
        Vector6 dsvMdsij = Vector6::Zero();
        double s_vM_2 = 0;
        double dAsdDL = 0;

        double ds_eq_2_dDL = 0;
        double ds_eq_2_dfM = 0;

        double dR2dDL = 0;
        double RA = 0;
        double s_eq_prime = sqrt(1.5*double_contract(stilde2, stilde2));
        double s_eq_2 = s_eq_prime;
        Matrix6x6 nnt = Matrix6x6::Zero();
        Matrix6x6 Aijkl = Matrix6x6::Zero();
        Vector6 nij2 = 1.5*stilde2;
        if (s_eq_2 > 1e-12) {
            nij2 /= s_eq_2;
        }
        else {
            nij2 *= 0;
        }

        Vector6 bij = Vector6::Zero();

        double R2 = 0;

        Vector6 dsij_prime_dDL = Vector6::Zero();
        double B = 1;
        double residual = 1e99;
        unsigned iter = 0;
        while (residual > 1e-15) {
            ++iter;
            // print_at_time("Newton loop iter " + std::to_string(iter), time[1], noel, npt);
            double fM2 = state.fM() + DfM;
            sigma_2 = sigma_t;

            double dDL = 0;
            double dDfM_stress = 0;
            double dDfM_strain = 0;

            B = 1 + 3*G*params.R2()*DfM/params.sy0A();
            // s_eq_2 = (s_eq_prime - 3*G*(DL + params.R1()*DfM))/B;
            I1_2 = sigma_t[0] + sigma_t[1] + sigma_t[2] - 3*K*params.dV()*DfM;
            s_eq_2 = (s_eq_prime - 3*G*(DL + params.R1()*DfM))/B;
            // Calculates f and the derivative df/dDL,
            if (plastic) {
                dsij_prime_dDL = Vector6::Zero();
                ds_eq_2_dDL = -3*G;
                double sy0 = params.sy0M()*(state.fM() + DfM) + params.sy0A()*(1 - (state.fM() + DfM));
                R2 = (state.R() + params.b()*params.Q()*DL)/(1 + params.b()*DL);
                double sy_2 = sy0 + R2;
                Vector6 sij_prime = sij_t;
                double back_stress_correction = 0;
                for (unsigned i = 0; i != params.back_stresses(); ++i) {
                    double theta = 1./(1 + params.gamma(i)*DL);
                    sij_prime -= theta*state.back_stress_vector(i);
                    back_stress_correction += theta*params.Cm(i)*DL;
                    dsij_prime_dDL += params.gamma(i)*theta*theta*state.back_stress_vector(i);
                    ds_eq_2_dDL -= theta*theta*params.Cm(i);
                }
                s_eq_prime = sqrt(1.5*double_contract(sij_prime, sij_prime));
                s_eq_2 = (s_eq_prime - 3*G*(DL + params.R1()*DfM) - back_stress_correction)/B;
                nij2 = 1.5*sij_prime/s_eq_prime;
                ds_eq_2_dDL += double_contract(nij2, dsij_prime_dDL);
                ds_eq_2_dDL /= B;
                dR2dDL = params.b()/(1 + params.b()*DL)*(params.Q() - R2);
                dfdDL = ds_eq_2_dDL - dR2dDL;
                f = s_eq_2 + params.a()*I1_2 - sy_2;
            }
            RA = params.R1() + params.R2()*s_eq_2/params.sy0A();
            ds_eq_2_dfM = -3*G*RA/B;
            dfdDfM = ds_eq_2_dfM - 3*params.a()*K*params.dV() - (params.sy0M() - params.sy0A());
            sigma_2 -= K*params.dV()*delta_ij*DfM;
            Vector6 dsijdDfM = -K*params.dV()*delta_ij;

            if (s_eq_prime > 1e-12) {
                sigma_2 -= 2*G*(DL + RA*DfM)*nij2;
                dsijdDfM -= 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2;
            }
            Vector6 dsijdDL = -2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2;
            // Calculating the von Mises stress at step 2
            s = deviator(sigma_2);
            double J2 = 0.5*double_contract(s, s);
            s_vM_2 = sqrt(3*J2);
            dsvMdsij = Vector6::Zero();
            if (J2 > 1e-12)
                dsvMdsij = 1.5*s/sqrt(3*J2);

            // h_strain and derivatives of h_strain
            if (plastic) {
                // print_at_time("Entering plastic section", time[1], noel, npt);
                Sigma = I1_2/s_vM_2;
                double DI1 = 3*K*(de[0] + de[1] + de[2] - DfM*params.dV());
                Vector6 Dsij = static_cast<Vector6>(2*G*(deviator(static_cast<Vector6>(de)) - (DL + RA*DfM)*nij2));
                double DvM = 1.5/s_vM_2*double_contract(deviator(sigma_2), Dsij);
                DSigma = Sigma*(DI1/I1_2 - DvM/s_vM_2);
                double n = params.n();
                double dSigmadDL = -Sigma/s_vM_2*double_contract(dsvMdsij, dsijdDL);
                double dSigmadDfM = 1/s_vM_2*(-3*K*params.dV() - Sigma*double_contract(dsvMdsij, dsijdDfM));
                double dDSigmadDL = -DSigma/DI1*double_contract(dsvMdsij, dsijdDL);
                double dDSigmadDfM = 1/DvM*(-3*K*params.dV() - DSigma*double_contract(dsvMdsij, dsijdDfM));
                fsb2 = 1 - (1 - state.fsb0())*exp(-params.alpha()*(state.ep_eff() + DL));
                dfsb2dDL = params.alpha()*(1 - state.fsb0())*exp(-params.alpha()*(state.ep_eff() + DL));

                double c = params.alpha()*params.beta()*n*(1 - fsb2)*pow(fsb2, n - 1);
                double dcdDL = n*params.alpha()*params.beta()*pow(fsb2, n - 2)*(n*(1 - fsb2) - 1)*dfsb2dDL;
                double Gamma = params.g0() - params.g1()*(temp-params.M_sigma())/(params.M_d() - params.M_sigma())
                               + params.g2()*Sigma;
                norm_drivning_force = (Gamma - params.g_mean())/params.g_std();
                P = 0.5*(1 + erf(norm_drivning_force));

                As = c*P;
                pdf = normal_pdf(norm_drivning_force)/params.g_std();
                Bs = params.g2()*params.beta()*pow(fsb2, n)*pdf*(DSigma > 0);
                dAsdDL = dcdDL*P + c*params.g2()*pdf*dSigmadDL;
                double dAsdfM = c*params.g2()*pdf*dSigmadDfM;

                double dBsdDL = (pdf*n*pow(fsb2, n - 1)*dfsb2dDL
                                 - Bs*norm_drivning_force/params.g_std()*params.g2()*dSigmadDL)*(DSigma > 0);
                double dBsdfM = -Bs*norm_drivning_force/params.g_std()*params.g2()*dSigmadDfM;
                h_strain = (1 - fM2)*(1 - state.other_phases())*(As*DL + Bs*DSigma) - DfM_strain;
                dh_straindDL =
                        (1 - fM2)*(1 - state.other_phases())*(As + DL*dAsdDL + Bs*dDSigmadDL + DSigma*dBsdDL);
                dh_strainDfM = -(As*DL + Bs*DSigma) +
                               (1 - fM2)*(1 - state.other_phases())*(DL*dAsdfM + Bs*dDSigmadDfM + DSigma*dBsdfM);
                // print_at_time("Plastic section done", time[1], noel, npt);
                // print_for_position("Plasticity done", "", noel, npt);
            }

            if (stress_transformations) {
                // print_for_position("Stress transformation: ", "", noel, npt);
                h_stress = stress_transformation_function(sigma_2, temp, params, state, fM2);
                bij = params.a1()*delta_ij;
                if (J2 > 1e-12) {
                    bij += params.a2()*dsvMdsij + params.a3()*(contract(s, s) - 2./3*J2*delta_ij);
                }
                double exp_fun = 1 - (h_stress + fM2);
                bij *= exp_fun*params.k();
                dh_stressDfM = double_contract(bij, dsijdDfM) - 1;

                // print_for_position("h_stress: ", h_stress, noel, npt);
                // print_for_position("dh_stressDfM: ", dh_stressDfM, noel, npt);
            }

            if (!plastic) {
                dDfM_stress = h_stress/dh_stressDfM;
            }
            else {
                if (!stress_transformations && !strain_transformations) {
                    // print_for_position("Only plastic", "", noel, npt);
                    // print_for_position("f: ", f, noel, npt);
                    // print_for_position("dfdDL: ", dfdDL, noel, npt);
                    if (dfdDL != 0) {
                        dDL = f/dfdDL;
                    }
                    else {
                        pnewdt = 0.25;
                        return;
                    }
                }
                else if (!stress_transformations) {
                    // print_for_position("Strain transformations", "", noel, npt);
                    double det = dfdDL*(dh_strainDfM - 1) - dfdDfM*dh_straindDL;
                    dDL = ((dh_strainDfM - 1)*f - dfdDfM*h_strain)/det;
                    dDfM_strain = (-dh_straindDL*f + dfdDL*h_strain)/det;
                }
                else if (!strain_transformations) {
                    double det = dfdDL*dh_stressDfM - dfdDfM*dh_stressDfM;
                    dDL = (dh_stressDfM*f - dfdDfM*h_stress)/det;
                    dDfM_stress = (-dh_stressDfM*f + dfdDL*h_stress)/det;
                }
                else {
                    double a = dfdDL;
                    double b = dfdDfM;
                    double c = dh_stressDL;
                    double d = dh_stressDfM;
                    double e = dh_straindDL;
                    double g = dh_strainDfM;
                    double det = a*d - c*b;   // Pseudo determinant appearing in the inverse
                    dDL = (d*f - b*h_stress)/det;
                    dDfM_stress = -(c + d*e - c*g)/det*f + (a + b*e - a*g)/det*h_stress + h_strain;
                    dDfM_strain = (d*e - c*g)/det*f - (b*e - a*g)/det*h_stress - h_strain;
                }
                if ((DL - dDL < 0 || DfM_strain - dDfM_strain < 0) && stress_transformations) {
                    DL = 0;
                    dDL = 0;
                    dDfM_strain = 0;
                    DfM_strain = 0;
                    dDfM_stress = h_stress/dh_stressDfM;

                }
                if (DfM_stress - dDfM_stress < 0) {
                    DfM_stress = 0;
                    dDfM_stress = 0;
                    double det = dfdDL*(dh_strainDfM - 1) - dfdDfM*dh_straindDL;
                    dDL = ((dh_strainDfM - 1)*f - dfdDfM*h_strain)/det;
                    dDfM_strain = (-dh_straindDL*f + dfdDL*h_strain)/det;
                }
            }
            // print_for_position("Done with updating iteration ", iter, noel, npt);

            DL -= dDL;
            DfM_stress -= dDfM_stress;
            // print_for_position("dDfM_stress", dDfM_stress, noel, npt);
            DfM_strain -= dDfM_strain;
            DfM = DfM_stress + DfM_strain;
            residual = abs(dDL) + abs(dDfM_stress) + abs(dDfM_strain);
            if (iter > 25) {
                pnewdt = 0.25;
                return;
            }
        }

        // Updating state variables
        // print_for_position("Increase in martensite", DfM, noel, npt);
        state.ep_eff() += DL;
        state.q_martensite() += DfM;
        state.austenite() -= DfM;
        state.fM_stress() += DfM_stress;
        state.fM_strain() += DfM_strain;
        state.R() = R2;
        state.fsb() = fsb2;
        stress_vec = sigma_2;

        if (params.kinematic_hardening()) {
            state.total_back_stress() = Vector6::Zero();
            for (unsigned i = 0; i != params.back_stresses(); ++i) {
                state.back_stress_vector(i) += 2./3*params.Cm(i)*DL*nij2;
                state.back_stress_vector(i) /= (1 + params.gamma(i)*DL);
                state.total_back_stress() += state.back_stress_vector(i);
            }
        }
        // print_at_time("Calculating tangent", time[1], noel, npt);
        nnt = nij2*nij2.transpose();
        Aijkl = J - 2./3*nnt;
        D_alg = Del;
        if (s_eq_prime > 1e-12) {
            D_alg -= 6*G*G*(DL + RA*DfM)/s_eq_prime*Aijkl;
        }

        double A = dR2dDL - ds_eq_2_dDL;

        Vector6 Lekl = Vector6::Zero();
        Vector6 Lskl = Vector6::Zero();
        Vector6 Fekl = Vector6::Zero();
        Vector6 Fskl = Vector6::Zero();


        if (DfM_stress > 0) {
            Fskl = bij;
            if (DL > 0) {
                Lekl = 1./A/B*(2*G*nij2 + params.a()*K*delta_ij);
                Lskl = 1./A*dfdDfM*bij;
            }
        }

        else {
            double B1 = (1 - state.fM())*(1 - state.other_phases())/(1 + As*DL + Bs*DSigma);
            double B2 = B1*(As + DL*dAsdDL);
            Lekl = (2*G*nij2 + params.a()*K*delta_ij)/(A*B - B*B2*dfdDfM);

            Vector6 Gkl =
                    B1*Bs*Sigma*(1 - norm_drivning_force/params.g_std()*params.g2())*(delta_ij/I1_2 - 1.5*s/s_vM_2);
            Lskl = Gkl*dfdDfM/(A + B2*dfdDfM);
            Fekl = B2*Lekl;
            Fskl = Gkl*(1 + dfdDfM/(A + B2*dfdDfM)*B2);
        }

        if (DL > 0) {
            D_alg -= 2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2*Lekl.transpose();
        }

        if (DfM > 0) {
            D_alg -= 4*G*G/B*DfM*params.R2()/params.sy0A()*nnt;
            Matrix6x6 Bijkl = I;
            if (DL > 0) {
                D_alg -= 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2*Fekl.transpose()
                         - K*params.dV()*delta_ij*Fekl.transpose();
                Bijkl += 2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2*Lskl.transpose();
            }
            Bijkl += 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2*Fskl.transpose()
                     + K*params.dV()*delta_ij*Fskl.transpose();

            D_alg = Bijkl.inverse()*D_alg;
        }
        if (std::fetestexcept(FE_INVALID)) {
            std::cout << "DL: " << DL << std::endl;
            std::cout << "nij: " << nij2.transpose().format(CleanFmt) << std::endl;
            std::cout << "s_eq_2: " << s_eq_2 << std::endl;
            std::cout << "s_eq_prime: " << s_eq_prime << std::endl;
            std::cout << "A: " << A << std::endl;
            std::cout << "B: " << B << std::endl;
            xit_();
        }
    }

}
