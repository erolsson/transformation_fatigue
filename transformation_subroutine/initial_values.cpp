//
// Created by erolsson on 06/11/2019.
//

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "transformation_subroutine_header.h"

struct HeatTreatmentDataCompare {
public:
    bool operator()(const HeatTreatmentData& a, const HeatTreatmentData& b) const {
        if (a.element == b.element) {
            return a.gauss_point < b.gauss_point;
        }
        return a.element < b.element;
    }
};

std::array<std::size_t, 8> gp_order_x_neg = {2, 1, 4, 3, 6, 5, 8, 7};

std::mutex part_info_mutex;

std::vector<HeatTreatmentData> heat_treatment_data;

std::string get_fortran_string(void (*fortran_func)(char*, int&, int)){
    char out_char[256];
    int out_len = 0;
    fortran_func(out_char, out_len, 256);
    std::string string = std::string(out_char, out_char+out_len);
    return string;

}

extern "C" void uexternaldb_(const int* lop, const int* lrestart, const double* time, const double* dtime,
        const int* kstep, const int* kinc) {
    if (*lop == 0) {
        // Beginning of analysis
        std::string run_directory = get_fortran_string(getoutdir_);
        std::string job_name = get_fortran_string(getjobname_);
        std::string data_file_name = job_name + ".htd";
        std::ifstream data_file(run_directory + "/" + data_file_name);
        std::string data_line;
        while (getline(data_file, data_line)) {
            std::vector<std::string> line_data;
            std::stringstream line(data_line);
            std::string val;
            while (getline(line, val, ',')) {
                line_data.push_back(val);
            }

            std::size_t stress_components = line_data.size() - 11;
            HeatTreatmentData data {};
            data.element = stoi(line_data[0]);
            data.gauss_point = stoi(line_data[1]);
            data.stress = Eigen::VectorXd(stress_components);
            for( unsigned i = 0; i != stress_components; ++i ) {
                data.stress(i) = stod(line_data[2+i]);
            }

            for (unsigned i = 0; i != 9; ++i) {
                data.phase_data.push_back(stod(line_data[2 + stress_components + i]));
            }
            heat_treatment_data.push_back(data);
        }
    }
}

std::size_t reorder_gauss_pt(std::size_t gp, std::string part_name) {

    if (part_name.find("x_neg") != std::string::npos) {
        return gp_order_x_neg[gp - 1];
    }
    return gp;
}

std::pair<std::size_t, std::string> user_model_data(int noel, const double* coords) {
    int user_elem_number = 0;
    {
        std::lock_guard<std::mutex> lock(part_info_mutex);
        getelemnumberuser_(noel, user_elem_number);
    }
    std::string part_name;
    if (coords[0] < 0) {
        part_name = "x_neg";
    }
    else {
        part_name = "x_pos";
    }
    return std::make_pair(user_elem_number, part_name);
}

std::vector<HeatTreatmentData>::iterator find_heat_treatment_data(int noel, int npt) {

    HeatTreatmentData tmp;
    tmp.element = noel;
    tmp.gauss_point = npt;
    auto it = std::lower_bound(heat_treatment_data.begin(), heat_treatment_data.end(), tmp,
                               HeatTreatmentDataCompare());
    if (tmp.element != it->element || tmp.gauss_point != it->gauss_point) {
        std::cout << "Element numbering in model does not correspond to the numbering in the data file, exiting!";
        std::cout << "Element number is " << tmp.element << " and Gauss point is " << tmp.gauss_point << std::endl;
        if (it == heat_treatment_data.end()) {
            std::cout << "The element is not found!" << std::endl;
        }
        else {
            std::cout << "Found element has element number " << it->element << " and Gauss point " << it->gauss_point
                      << std::endl;
        }
        std::abort();
    }
    return it;
}

extern "C" void sdvini_(double* statev, const double* coords, const int& nstatev, const int& ncrds, const int& noel,
                        const int& npt, const int& layer, const int& kspt) {
    auto user_data = user_model_data(noel, coords);
    std::size_t gp = reorder_gauss_pt(npt, user_data.second);
    auto it = find_heat_treatment_data(user_data.first, gp);
    statev[0] = 0.;
    for (unsigned i = 0; i != 9; ++i) {
        if (i != 3) {
            statev[i+1] = it->phase_data[i];
        }
        else {
            double HRC = it->phase_data[i];
            double HV = (223*HRC - 14500)/(100-HRC);
            statev[i+1] = HV;
        }
    }
    for( unsigned i = 10; i != nstatev; ++i) {
        statev[i] = 0;
    }
}

extern "C" void sigini_(double* sigma, const double* coords, const int& ntens, const int& ncords, const int& noel,
                        const int& npt, const int& layer, const int& kspt, const int& rebar, const char* names) {
    auto user_data = user_model_data(noel, coords);
    std::size_t gp = reorder_gauss_pt(npt, user_data.second);
    auto it = find_heat_treatment_data(user_data.first, gp);
    for (unsigned i = 0; i != ntens; ++i) {
        if (user_data.second.find("x_neg") != std::string::npos && (i == 3 || i == 4)) {
            sigma[i] = -(it->stress(i));
            std::cout << "Neg stress" << std::endl;
        }
        else {
            sigma[i] = it->stress(i);
            std::cout << "Pos stress" << std::endl;
        }
    }
}