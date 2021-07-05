//
// Created by erolsson on 19/12/2019.
//


#ifndef TRANSFORMATION_FATIGUE_TRANSFORMATION_SUBROUTINE_HEADER_H
#define TRANSFORMATION_FATIGUE_TRANSFORMATION_SUBROUTINE_HEADER_H

#include <array>
#include <vector>
#include <Eigen/Dense>


extern "C" void getoutdir_(char* outdir, int&, int);
extern "C" void getjobname_(char* outdir, int&, int);
extern "C" void getpartinfoc_(char* name, int& name_len, const int& num, const int& jtype, int& user_num, int& error);
extern "C" void getelemnumberuser_(const int& num, int& user_num);
extern "C" void getpartname_(char* outdir, int&, int);
extern "C" void rotsig_(double* S, double* R, double* Srot, const int& lstr, const int& ndi,
                        const int& nshr);
extern "C" void xit_();


#endif //TRANSFORMATION_FATIGUE_TRANSFORMATION_SUBROUTINE_HEADER_H

