#ifndef PYBIND11_H
#define PYBIND11_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
#endif

#ifndef EUCHAR_LIBUTILS_H
#define EUCHAR_LIBUTILS_H

using namespace std;

int sum_bool_2d(vector<vector<bool>> matrix);

vector<vector<int>> pad_2d(const vector<vector<int>> &image, int M);

vector<vector<bool>> threshold_image_2d(const vector<vector<int>> &image, int value);

vector<vector<bool>> elementwise_AND_2d(const vector<vector<bool>> &image1, const vector<vector<bool>> &image2);

vector<vector<int>> neigh_pixel_2d(const vector<vector<int>> &padded, size_t i , size_t j);

vector<vector<bool>> binary_neigh_pixel_2d(const vector<vector<int>> &padded, size_t i , size_t j , int pixel_value);

int char_binary_image_2d(vector<vector<bool>> input);

vector<vector<bool>> neigh_2d_from_number(size_t num);

size_t number_from_neigh_2d(vector<vector<bool>> neigh);

vector<short> vector_of_euler_changes_2d();

int sum_bool_3d(vector<vector<vector<bool>>> matrix);

vector<vector<vector<int>>> pad_3d(py::array_t<int> input, int M);

vector<vector<vector<bool>>> threshold_image_3d(const vector<vector<vector<int>>> &image, int value);

vector<vector<vector<bool>>> elementwise_AND_3d(const vector<vector<vector<bool>>> &image1, const vector<vector<vector<bool>>> &image2);

vector<vector<vector<int>>> neigh_voxel_3d(const vector<vector<vector<int>>> &padded, size_t i , size_t j , size_t k);

vector<vector<vector<bool>>> binary_neigh_voxel_3d(const vector<vector<vector<int>>> &padded, size_t i , size_t j , size_t k, int voxel_value);

int char_binary_image_3d(vector<vector<vector<bool>>> input);

vector<vector<vector<bool>>> neigh_3d_from_number(size_t num);

size_t number_from_neigh_3d(vector<vector<vector<bool>>> neigh);

vector<short> vector_of_euler_changes_3d(size_t offset=0, size_t max_value=67108864);

size_t dim_simplex(const vector<int> &simplex);

vector<bool> filter_parametrization_on(const vector<int> &vertices, const vector<vector<int>> &simplices, const vector<double> &param);

#endif