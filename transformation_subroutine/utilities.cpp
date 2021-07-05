//
// Created by erolsson on 18/11/2020.
//

#include <iostream>
#include <mutex>
#include <string>

#include "utilities.h"

std::mutex part_info_mutex;

std::string get_fortran_string(void (*fortran_func)(char*, int&, int)){
    char out_char[256];
    int out_len = 0;
    fortran_func(out_char, out_len, 256);
    return std::string(out_char, out_char+out_len);
}

std::string get_part_name(int element_number) {
    int name_len = 80;
    char name[name_len];
    {
        std::lock_guard <std::mutex> lock(part_info_mutex);
        getpartname_(element_number, name, name_len);
    }
    std::string data(name);
    auto first_space = data.find(" ");
    return std::string(data.begin(), data.begin() + first_space);
}

int get_user_element_number(int element_number) {
    int user_number = 1;
    {
        std::lock_guard <std::mutex> lock(part_info_mutex);
        getelemnumberuser_(element_number, user_number);
    }
    return user_number;
}

std::string get_job_name() {
    int len = 255;
    int job_name_len = 0;
    char name[len];
    {
        std::lock_guard <std::mutex> lock(part_info_mutex);
        getjobname_(name, job_name_len, len);
    }
    return std::string(name, name + job_name_len);
}

std::string get_job_directory() {
    int len = 255;
    int directory_name_len = 0;
    char directory[len];
    {
        std::lock_guard <std::mutex> lock(part_info_mutex);
        getoutdir_(directory, directory_name_len, len);
    }
    return std::string(directory, directory + directory_name_len);
}
