#ifndef PYBIND11_H
#define PYBIND11_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
#endif

#ifndef EUCHAR_LIBCURVE_H
#define EUCHAR_LIBCURVE_H

using namespace std;

vector<int> naive_image_2d(vector<vector<int>> image, int M);

vector<int> image_2d(const vector<vector<int>> & image, const vector<int> &vector_euler_changes, int M);

vector<int> naive_image_3d(vector<vector<vector<int>>> image, int M);

vector<int> image_3d(py::array_t<int> input, const vector<int> &vector_euler_changes, int M);

vector<int> filtration_2d(const vector<vector<int>> &simplices, const vector<double> &param, vector<double> &bins);

#endif