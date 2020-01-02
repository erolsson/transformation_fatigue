//
// Created by erolsson on 19/12/2019.
//


#ifndef TRANSFORMATION_FATIGUE_TRANSFORMATION_SUBROUTINE_HEADER_H
#define TRANSFORMATION_FATIGUE_TRANSFORMATION_SUBROUTINE_HEADER_H

#include <array>
#include <vector>
#include <Eigen/Dense>

struct HeatTreatmentData {
    std::size_t element = 0;
    std::size_t gauss_point = 0;
    Eigen::VectorXd stress;
    std::vector<double> phase_data;
};

extern "C" void getoutdir_(char* outdir, int&, int);
extern "C" void getjobname_(char* outdir, int&, int);
extern "C" void getpartinfoc_(char* name, const int& num, const int& jtype, int& user_num, int& error);
extern "C" void getelemnumberuser_(const int& num, int& user_num);
extern "C" void getpartname_(char* outdir, int&, int);

#endif //CASE_HARDENING_TOOLBOX_COOLING_SUBROUTINE_HEADER_H

