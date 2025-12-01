#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <bits/stdc++.h>
#include "M_function.h"

surface M_init(const std::string discharging_filename, const double nominal_capacity, const double R_id, std::vector<double>& a_1, std::vector<double>& a_1_I, std::vector<double>& a_2, std::vector<double>& a_2_I, int S, int P);
void exportcurveToCSV(const std::string& filename, std::vector<double> x, std::vector<double> y);
std::vector<std::vector<double>> readCSV(const std::string& filename);

