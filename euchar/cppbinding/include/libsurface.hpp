#ifndef PYBIND11_H
#define PYBIND11_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
#endif

#ifndef EUCHAR_LIBSURFACE_H
#define EUCHAR_LIBSURFACE_H

using namespace std;

vector<vector<int>> naive_images_2d(const vector<vector<int>> &image1, const vector<vector<int>> &image2, int M1, int M2);

vector<vector<int>> images_2d(const vector<vector<int>> &image1, const vector<vector<int>> &image2, const vector<int> &vector_euler_changes, int M1, int M2);

vector<vector<int>> naive_images_3d(const vector<vector<vector<int>>> &image1, const vector<vector<vector<int>>> &image2, int M1, int M2);

vector<vector<int>> images_3d(py::array_t<int> input1, py::array_t<int> input2, vector<int> vector_euler_changes, int M1 = 255, int M2 = 255);

vector<vector<int>> bifiltration(vector<int>    dim_simplices, vector<double> parametrization1, vector<double> parametrization2, vector<double> bins1, vector<double> bins2);

vector<vector<int>> bifiltration_2d(vector<double> unsorted_par_vertices, vector<double> sorted_par_simplices, vector<vector<int>> sorted_simplices, size_t nbins1, size_t nbins2, vector<double> minmax1, vector<double> minmax2);

vector<vector<int>> bifiltration_3d(vector<double> unsorted_par_vertices, vector<double> sorted_par_simplices, vector<vector<int>> sorted_simplices, size_t nbins1, size_t nbins2, vector<double> minmax1, vector<double> minmax2);

#endif