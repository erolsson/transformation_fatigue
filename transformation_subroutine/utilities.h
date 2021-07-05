//
// Created by erolsson on 18/11/2020.
//

#ifndef CASE_HARDENING_TOOLBOX_UTILITIES_H
#define CASE_HARDENING_TOOLBOX_UTILITIES_H

#include <string>
#include <vector>

#include <Eigen/Dense>

// Small struct for handling heat treatment simulation results in each gauss point
struct HeatTreatmentData {
    std::size_t element = 0;
    std::size_t gauss_point = 0;
    Eigen::VectorXd stress;
    std::vector<double> phase_data;
};

struct HeatTreatmentDataCompare {
public:
    bool operator()(const HeatTreatmentData& a, const HeatTreatmentData& b) const {
        if (a.element == b.element) {
            return a.gauss_point < b.gauss_point;
        }
        return a.element < b.element;
    }
};

inline
Eigen::VectorXd rotate_stress_state_z(const Eigen::VectorXd& stress_state, double angle) {
    Eigen::Matrix3d S;
    S << stress_state(0), stress_state(3), stress_state(4),
         stress_state(3), stress_state(1), stress_state(5),
         stress_state(4), stress_state(5), stress_state(2);
    Eigen::Matrix3d Q;
    Q << cos(angle), sin(angle), 0,
        -sin(angle), cos(angle), 0,
        0, 0, 1;
    S = Q*S*Q.transpose();
    Eigen::VectorXd s = Eigen::VectorXd::Zero(6);
    s(0) = S(0,0);
    s(1) = S(1,1);
    s(2) = S(2,2);
    s(3) = S(0,1);
    s(4) = S(0,2);
    s(5) = S(1,2);
    return s;
}

// Functions for a smoother interface with abaqus fortran interface
std::string get_part_name(int element_number);
int get_user_element_number(int element_number);
std::string get_job_name();
std::string get_job_directory();

// External declarations of abaqus fortran functions used in the above functions

extern "C" void getoutdir_(char* outdir, int&, int);
extern "C" void getjobname_(char* outdir, int&, int);
extern "C" void getpartinfoc_(char* name, int& name_len, const int& num, const int& jtype, int& user_num, int& error);
extern "C" void getelemnumberuser_(const int& num, int& user_num);
extern "C" void getpartname_(int& num, char* name, int& name_len);

#endif //CASE_HARDENING_TOOLBOX_UTILITIES_H
