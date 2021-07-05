//
// Created by erolsson on 02/07/2021.
//

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>

#include "utilities.h"

std::map<std::string, std::vector<HeatTreatmentData>> heat_treatment_data {
        {"coarse", { }},
        {"dense", { }}
};

extern "C" void uexternaldb_(const int* lop, const int* lrestart, const double* time, const double* dtime,
                             const int* kstep, const int* kinc) {
    if (*lop == 0) {
        // Beginning of analysis
        std::string run_directory = get_job_directory();
        std::string job_name = get_job_name();
        auto cd_pos = job_name.find("cd=");
        std::string case_depth = std::string(job_name.begin() + cd_pos + 3, job_name.begin() + cd_pos + 5);
        std::array<std::string, 2> teeth = {"coarse", "dense"};
        for(const auto tooth: teeth) {
            std::string heat_treatment_filename = "cd=" + case_depth + "_" + tooth + ".htd";
            std::ifstream data_file(run_directory + "/" + heat_treatment_filename);
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
                heat_treatment_data[tooth].push_back(data);

            }
        }
    }
}


std::vector<HeatTreatmentData>::iterator find_heat_treatment_data(int user_element_number, int gauss_point,
                                                                  const std::string& part_type) {

    HeatTreatmentData tmp;
    tmp.element = user_element_number;
    tmp.gauss_point = gauss_point;
    auto it = std::lower_bound(heat_treatment_data[part_type].begin(), heat_treatment_data[part_type].end(), tmp,
                               HeatTreatmentDataCompare());
    if (tmp.element != it->element || tmp.gauss_point != it->gauss_point) {
        std::cout << "Element numbering in model does not correspond to the numbering in the data file, exiting!"
                  << std::endl;
        std::cout << "Element number is " << tmp.element << " and Gauss point is " << tmp.gauss_point << std::endl;
        if (it == heat_treatment_data[part_type].end()) {
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


std::size_t reorder_gauss_pt(std::size_t gp, const std::string& sign) {
    std::array<std::size_t, 8> gp_order_x_neg = {2, 1, 4, 3,
                                                 6, 5, 8, 7};
    if (sign == "neg") {
        return gp_order_x_neg[gp - 1];
    }
    return gp;
}

extern "C" void sigini_(double* sigma, const double* coords, const int& ntens, const int& ncords, const int& noel,
                        const int& npt, const int& layer, const int& kspt, const int& rebar, const char* names) {
    Eigen::VectorXd stress_vec = Eigen::VectorXd::Zero(ntens);
    std::string part_name = get_part_name(noel);
    std::transform(part_name.begin(), part_name.end(), part_name.begin(),
                   [](unsigned char c){ return std::tolower(c);});
    std::string part_type = "jaw";
    auto part_type_pos = part_name.find("dense");
    if (part_type_pos != std::string::npos) {
        part_type = "dense";
    }
    else {
        part_type_pos = part_name.find("coarse");
        if (part_type_pos != std::string::npos) {
            part_type = "coarse";
        }
    }
     if (part_type == "coarse" || part_type == "dense") {
        int user_element_number = get_user_element_number(noel);
        std::string sign = std::string(part_name.begin() + (part_type + "_tooth_").size(),
                           part_name.begin() + (part_type + "_tooth_").size() + 3);
        std::size_t gp = reorder_gauss_pt(npt, sign);
        auto it = find_heat_treatment_data(user_element_number, gp, part_type);
        stress_vec = it->stress;
        if (sign == "neg") {
            stress_vec(3) *= -1;
            stress_vec(4) *= -1;
        }
        double angle = atan2(coords[1], coords[0]);
        stress_vec = rotate_stress_state_z(stress_vec, angle);
    }

    for (unsigned i = 0; i != ntens; ++i) {
        sigma[i] = stress_vec(i);
    }
}

extern "C" void sdvini_(double* statev, const double* coords, const int& nstatev, const int& ncrds, const int& noel,
                        const int& npt, const int& layer, const int& kspt) {
    std::string part_name = get_part_name(noel);
    std::transform(part_name.begin(), part_name.end(), part_name.begin(),
                   [](unsigned char c){ return std::tolower(c);});
    std::string part_type = "jaw";
    auto part_type_pos = part_name.find("dense");
    if (part_type_pos != std::string::npos) {
        part_type = "dense";
    }
    else {
        part_type_pos = part_name.find("coarse");
        if (part_type_pos != std::string::npos) {
            part_type = "coarse";
        }
    }
    if (part_type == "coarse" || part_type == "dense") {
        int user_element_number = get_user_element_number(noel);
        std::string sign = std::string(part_name.begin() + (part_type + "_tooth_").size(),
                                       part_name.begin() + (part_type + "_tooth_").size() + 3);
        std::size_t gp = reorder_gauss_pt(npt, sign);
        auto it = find_heat_treatment_data(user_element_number, gp, part_type);
        statev[0] = 0.;
        for (unsigned i = 0; i != 9; ++i) {
            statev[i+1] = it->phase_data[i];
        }
        for( unsigned i = 10; i != nstatev; ++i) {
            if (i == 12 || i == 13) {
                statev[i] = 0.23;       // Initial shear band fraction
            }
            else {
                statev[i] = 0;
            }
        }
    }
    else {
        for (unsigned i = 0; i != nstatev; ++i) {
            statev[i] = 0;
        }
    }
}