#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>     // py::array_t<>
#include <vector>
#include <algorithm>    // std::lower_bound,

#include "libutils.hpp"

#include <cstdio>

namespace py = pybind11;
using namespace std;

//================================================

vector<vector<int>> naive_images_2d(const vector<vector<int>> &image1,
                                    const vector<vector<int>> &image2,
                                    int M1, int M2)
{
    if ( image1.size() != image2.size() || image1[0].size() != image2[0].size())
        throw runtime_error("Images dimensions must be equal.");

    vector<vector<int>> surface(M1+1, vector<int>(M2+1, 0));

    // Loop over all pixel in original image
    for (int s = 0; s < M1+1; ++s) {
        for (int t = 0; t < M2+1; ++t) {
            vector<vector<bool>> th1 = threshold_image_2d(image1, s);
            vector<vector<bool>> th2 = threshold_image_2d(image2, t);
            vector<vector<bool>> K_st = elementwise_AND_2d(th1, th2);
            surface[s][t] = char_binary_image_2d(K_st);
        }
    }

    return surface;
}


//================================================

vector<vector<int>> images_2d(const vector<vector<int>> &image1,
                              const vector<vector<int>> &image2,
                              const vector<int> &vector_euler_changes,
                              int M1, int M2)
{
    if ( image1.size() != image2.size() || image1[0].size() != image2[0].size() )
        throw runtime_error("Images dimensions must be equal.");
    
    size_t numI = image1.size();
    size_t numJ = image1[0].size();
    
    vector<vector<int>> surface(M1+1, vector<int>(M2+1, 0));

    // Padded images
    vector<vector<int>> padded1 = pad_2d(image1, M1);
    vector<vector<int>> padded2 = pad_2d(image2, M2);

    // Loop over all pixel in original image
    for (size_t i = 1; i < numI+1; ++i) {
        for (size_t j = 1; j < numJ+1; ++j) {
            int pixel_value1 = padded1[i][j];
            vector<vector<int>> neigh1 = neigh_pixel_2d(padded1, i, j);
            vector<vector<bool>> binary_neigh1 = binary_neigh_pixel_2d(padded1, i, j, pixel_value1);

            vector<vector<int>> neigh2 = neigh_pixel_2d(padded2, i, j);
            
            /* Sort and make unique pixels in neigh of padded2.
            These can be used to obtain thresholded neigh of
            padded2 used to compute all the euler changes along a 
            row of the euler surface.*/
            vector<int> indices_row(10, M2+1);
            int central_pixel = neigh2[1][1];
            for (size_t p = 0; p < 3; ++p) {
                for (size_t q = 0; q < 3; ++q) {
                    int pixel_value2 = neigh2[p][q];
                    if (central_pixel <= pixel_value2)
                        indices_row[q + p*3] = pixel_value2;
                }
            }
            
            // Sort and make unique indices_row
            sort(indices_row.begin(), indices_row.end());
            auto last = unique(indices_row.begin(), indices_row.end());
            indices_row.erase(last, indices_row.end());

            for (size_t z = 0; z < indices_row.size()-1; ++z) {
                size_t start = static_cast<size_t>(indices_row[z]);
                size_t end   = static_cast<size_t>(indices_row[z+1]);
                
                vector<vector<bool>> binary_neigh2 = threshold_image_2d(neigh2, static_cast<int>(indices_row[z]));
                vector<vector<bool>> binary_neigh = elementwise_AND_2d(binary_neigh1, binary_neigh2);
                size_t num = number_from_neigh_2d(binary_neigh);
                int change = static_cast<int>(vector_euler_changes[num]);

                for (size_t ind = start; ind < end; ++ind) {
                    surface[static_cast<size_t>(pixel_value1)][ind] += change;
                }
                
            }    
        }
    }
    
    //Cumulative sum along columns
    vector<int> row_before = surface[0];
    for (size_t s = 1; s < M1+1; ++s) {
        for (size_t t = 0; t < M2+1; ++t) {
             surface[s][t] += row_before[t];
        }
        row_before = surface[s];
    }

    return surface;
}

//================================================

vector<vector<int>> naive_images_3d(const vector<vector<vector<int>>> &image1,
                                    const vector<vector<vector<int>>> &image2,
                                    int M1, int M2)
{
    if ( image1.size() != image2.size() || image1[0].size() != image2[0].size() || image1[0][0].size() != image2[0][0].size())
        throw runtime_error("Images dimensions must be equal.");

    vector<vector<int>> surface(M1+1, vector<int>(M2+1, 0));

    // Loop over all pixel in original image
    for (int s = 0; s < M1+1; ++s) {
        for (int t = 0; t < M2+1; ++t) {
            vector<vector<vector<bool>>> th1 = threshold_image_3d(image1, s);
            vector<vector<vector<bool>>> th2 = threshold_image_3d(image2, t);
            vector<vector<vector<bool>>> K_st = elementwise_AND_3d(th1, th2);
            surface[s][t] = char_binary_image_3d(K_st);
        }
    }

    return surface;
}

//================================================

vector<vector<int>> images_3d(py::array_t<int> input1,
                              py::array_t<int> input2,
                              vector<int> vector_euler_changes,
                              int M1 = 255, int M2 = 255)
{
    auto image1 = input1.unchecked<3>();
    auto image2 = input2.unchecked<3>();
    
    if ( image1.shape(0) != image2.shape(0) || image1.shape(2) != image2.shape(2) || image1.shape(2) != image2.shape(2))
        throw runtime_error("Images dimensions must be equal.");

    size_t numI = image1.shape(0);
    size_t numJ = image1.shape(1);
    size_t numK = image1.shape(2);
    
    vector<vector<int>> surface(M1+1, vector<int>(M2+1, 0));
    
    // Padded images
    vector<vector<vector<int>>> padded1 = pad_3d(input1, M1);
    vector<vector<vector<int>>> padded2 = pad_3d(input2, M2);

    // Loop over all pixel in original image
    for (size_t i = 1; i < numI+1; ++i) {
        for (size_t j = 1; j < numJ+1; ++j) {
            for (size_t k = 1; k < numK+1; ++k) {
                int voxel_value1 = padded1[i][j][k];
                vector<vector<vector<int>>> neigh1 = neigh_voxel_3d(padded1, i, j, k);
                vector<vector<vector<bool>>> binary_neigh1 = binary_neigh_voxel_3d(padded1, i, j, k, voxel_value1);

                vector<vector<vector<int>>> neigh2 = neigh_voxel_3d(padded2, i, j, k);

                /* Sort and make unique voxels in neigh of padded2.
                   These can be used to obtain thresholded neigh of
                   padded2 used to compute all the euler changes
                   along a row of the euler surface.*/ 
                vector<int> indices_row(28, M2+1);  // 28=3^3+1
                int central_voxel2 = padded2[i][j][k];
                
                for (size_t p = 0; p < 3; ++p) {
                    for (size_t q = 0; q < 3; ++q) {
                        for (size_t r = 0; r < 3; ++r) {
                            int voxel_value2 = neigh2[p][q][r];
                            if (central_voxel2 <= voxel_value2)
                                indices_row[r + q*3 + p*9] = voxel_value2;
                        }
                    }
                }
                
                // Sort and make unique indices_row
                sort(indices_row.begin(), indices_row.end());
                auto last = unique(indices_row.begin(), indices_row.end());
                indices_row.erase(last, indices_row.end());

                for (size_t z = 0; z < indices_row.size()-1; ++z) {
                    size_t start = static_cast<size_t>(indices_row[z]);
                    size_t end   = static_cast<size_t>(indices_row[z+1]);
                    
                    vector<vector<vector<bool>>> binary_neigh2 = threshold_image_3d(neigh2, static_cast<int>(indices_row[z]));
                    vector<vector<vector<bool>>> binary_neigh = elementwise_AND_3d(binary_neigh1, binary_neigh2);
                    size_t num = number_from_neigh_3d(binary_neigh);
                    int change = static_cast<int>(vector_euler_changes[num]);
                    for (size_t ind = start; ind < end; ++ind) {
                        surface[static_cast<size_t>(voxel_value1)][ind] += change;
                    }
                }
            }
        }
    }
    
    //Cumulative sum along columns
    vector<int> row_before = surface[0];
    for (size_t s = 1; s < M1+1; ++s) {
        for (size_t t = 0; t < M2+1; ++t) {
             surface[s][t] += row_before[t];
        }
        row_before = surface[s];
    }

    return surface;
}

//================================================

vector<vector<int>> bifiltration(vector<int>    dim_simplices,
                                 vector<double> parametrization1,
                                 vector<double> parametrization2,
                                 vector<double> bins1,
                                 vector<double> bins2)
{
    size_t num_rows = bins1.size();
    size_t num_cols = bins2.size();

    vector<vector<int>> euler_char_surface(num_rows, vector<int>(num_cols, 0));

    // we have 4 possible changes, from dim 0 to dim 3 maximum
    vector<int> possible_changes = {1, -1, 1, -1};
    
    size_t number_simplices = dim_simplices.size();

    
    // loop on simplices
    for (size_t i = 0; i < number_simplices; ++i) {
        size_t dim_simplex = dim_simplices[i];
        double par1 = parametrization1[i];
        double par2 = parametrization2[i];
        
        // find bins indexes
        vector<double>::iterator lower1;
        lower1 = lower_bound(bins1.begin(), bins1.end(), par1);
        size_t index_simplex_in_bins1 = lower1 - bins1.begin();

        vector<double>::iterator lower2;
        lower2 = lower_bound(bins2.begin(), bins2.end(), par2);
        size_t index_simplex_in_bins2 = lower2 - bins2.begin();
        
        // update euler changes along rows
        for (size_t incr = 0; incr < num_cols - index_simplex_in_bins2; ++incr) {
            euler_char_surface[index_simplex_in_bins1][index_simplex_in_bins2+incr] += possible_changes[dim_simplex];
        }
    }

    //Cumulative sum along columns
    vector<int> row_before = euler_char_surface[0];
    for (size_t s = 1; s < num_rows; ++s) {
        for (size_t t = 0; t < num_cols; ++t) {
             euler_char_surface[s][t] += row_before[t];
        }
        row_before = euler_char_surface[s];
    }

    return euler_char_surface;
}

