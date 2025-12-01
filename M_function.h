#include"pch.h"
#include <fstream>
#include <iostream>
#include <cmath>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef K::Point_2 Point;

typedef struct surface_s{
  std::map<Point, double> function_values;
  Delaunay dt;
} surface;

surface s_init(std::vector<double> &all_currents, std::vector<double> &all_bks, std::vector<double> &all_vs);
double interpolate(double x, double y,const surface& m);
void exportSurfaceToCSV(const std::string& filename, 
                        double x_min, double x_max, int x_steps,
                        double y_min, double y_max, int y_steps,const surface& M);

