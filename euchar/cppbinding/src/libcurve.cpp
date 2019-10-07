#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>     // py::array_t<>
#include <numeric>      // std::partial_su, std::iota
#include <algorithm> // std::copy, std::lower_bound
#include <vector>
#include "libutils.hpp"
#include <cstdio>

namespace py = pybind11;

using namespace std;

//================================================

vector<int> naive_image_2d(vector<vector<int>> image, int M)
{
    vector<int> ecc(M+1, 0);

    for (size_t i = 0; i < M+1; ++i) {
        vector<vector<bool>> thresh_image = threshold_image_2d(image, static_cast<int>(i));
        ecc[i] = char_binary_image_2d(thresh_image);
    }

    return ecc;
}

//================================================

vector<int> image_2d(const vector<vector<int>> & image,
                     const vector<int> &vector_euler_changes,
                     int M)
{
    size_t numI = image.size();
    size_t numJ = image[0].size();

    vector<int> ecc(M+1, 0);

    vector<vector<int>> padded = pad_2d(image, M);

    for (size_t i = 1; i < numI+1; ++i) {
        for (size_t j = 1; j < numJ+1; ++j) {
            int pixel_value = padded[i][j];
            vector<vector<bool>> binary_neigh_3_3 = binary_neigh_pixel_2d(padded, i, j, pixel_value);
            size_t num = number_from_neigh_2d(binary_neigh_3_3);
            
            ecc[static_cast<size_t>(pixel_value)] += static_cast<int>(vector_euler_changes[num]);
        }
    }

    // Cumulative sum of euler changes
    partial_sum(ecc.begin(), ecc.end(), ecc.begin());

    return ecc;
}


//================================================

vector<int> naive_image_3d(vector<vector<vector<int>>> image, int M)
{
    vector<int> ecc(M+1, 0);

    for (size_t i = 0; i < M+1; ++i) {
        vector<vector<vector<bool>>> binary_thresh = threshold_image_3d(image, static_cast<int>(i));
        ecc[i] = char_binary_image_3d(binary_thresh);
    }

    return ecc;
}


//================================================

vector<int> image_3d(py::array_t<int> input,
                     const vector<int> &vector_euler_changes,
                     int M)
{
    auto image = input.unchecked<3>();
    
    size_t numI = image.shape(0);
    size_t numJ = image.shape(1);
    size_t numK = image.shape(2);

    vector<int> ecc(M+1, 0);

    vector<vector<vector<int>>> padded = pad_3d(input, M);

    for (size_t i = 1; i < numI+1; ++i) {
        for (size_t j = 1; j < numJ+1; ++j) {
            for (size_t k = 1; k < numK+1; ++k) {
                int voxel_value = padded[i][j][k];
                vector<vector<vector<bool>>> binary_neigh_3_3_3 = binary_neigh_voxel_3d(padded, i, j, k, voxel_value);
                size_t num = number_from_neigh_3d(binary_neigh_3_3_3);
            
            ecc[static_cast<size_t>(voxel_value)] += static_cast<int>(vector_euler_changes[num]);
            }    
        }
    }

    // Cumulative sum of euler changes
    partial_sum(ecc.begin(), ecc.end(), ecc.begin());

    return ecc;
}

//================================================

vector<int> filtration_2d(const vector<vector<int>> &simplices,
                          const vector<double> &param,
                          vector<double> &bins)
{
    size_t nbins(bins.size());
    vector<int> changes(nbins, 0);
    vector<int> ecc(nbins-1, 0);
    
    if (simplices.size() != param.size()) {
        cout << "Simplices and the parametrization must have";
        cout << " same number of elements." << endl;
        return ecc;
    }

    // possible changes due to addition of vertex edge or triangle
    vector<int> possible_changes{1, -1, 1};

    // loop on simplices and update euler curve
    for (size_t k=0; k < simplices.size(); k++) {
        vector<int> simpk = simplices[k];
        size_t dim = dim_simplex(simpk);

        vector<double>::iterator lower;
        lower = lower_bound(bins.begin(), bins.end(), param[k]);
        changes[(lower - bins.begin())] += possible_changes[dim];

    }

    int c = changes[0];
    for (size_t index = 0; index < ecc.size(); index++) {
        ecc[index] = c + changes[index+1];
        c = ecc[index];
    }
    
    return ecc;
}

