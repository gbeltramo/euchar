#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>     // py::array_t<>
#include <vector>
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
        printf("i = %zu | ", i);
        for (size_t j = 1; j < numJ+1; ++j) {
            printf("j = %zu | ", j);
            for (size_t k = 1; k < numK+1; ++k) {
                int voxel_value1 = padded1[i][j][k];
                vector<vector<vector<int>>> neigh1 = neigh_voxel_3d(padded1, i, j, k);
                vector<vector<vector<bool>>> binary_neigh1 = binary_neigh_voxel_3d(padded1, i, j, k, voxel_value1);

                vector<vector<vector<int>>> neigh2 = neigh_voxel_3d(padded2, i, j, k);

                /* Sort and make unique voxels in neigh of padded2.
                   These can be used to obtain thresholded neigh of
                   padded2 used to compute all the euler changes
                   along a row of the euler surface.*/ 
                vector<int> indices_row(10, M2+1);
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

vector<vector<int>> bifiltration_2d(vector<double> unsorted_par_vertices,
                                    vector<double> sorted_par_simplices,
                                    vector<vector<int>> sorted_simplices,
                                    size_t nbins1, size_t nbins2,
                                    vector<double> minmax1,
                                    vector<double> minmax2)
{
    // sorted indices of vertices
    vector<size_t> argvert(unsorted_par_vertices.size());
    iota(argvert.begin(), argvert.end(), 0);

    // map of vertices parametrization
    map<size_t, double> map_par_vertices;
    for (size_t k = 0; k < unsorted_par_vertices.size(); ++k)
        map_par_vertices[k] = unsorted_par_vertices[k];

    // sort indexes based on comparing values
    // in unsorted_par_vertices
    sort(argvert.begin(), argvert.end(),
        [&unsorted_par_vertices](size_t i1, size_t i2)
            {return unsorted_par_vertices[i1] < unsorted_par_vertices[i2];});

    // make bins1-->(vertices) and bins2-->(simplices)

    double step1((minmax1[1] - minmax1[0]) / nbins1);
    double step2((minmax2[1] - minmax2[0]) / nbins2);

    vector<double> bins1(nbins1+1);   // the last element of these
    vector<double> bins2(nbins2+1);   // is different from np.linspace()

    for (size_t row = 0; row < nbins1+1; row++)
        bins1[row] = step1 * row + minmax1[0];

    for (size_t col = 0; col < nbins2+1; col++)
        bins2[col] = step2 * col + minmax2[0];

    // map of simplices representing bi-filtration
    map<size_t, vector<size_t>> map_simplices;
    for (size_t index_simpl = 0; index_simpl < sorted_simplices.size(); index_simpl++) {
        vector<int> simpl(sorted_simplices[index_simpl]);

        double par_vertex_of_simpl(map_par_vertices[simpl[0]]);
        for (const int &vertex: simpl) {
            if (vertex == -1)
                continue;  // if simp is "sorted" use break

            if (map_par_vertices[vertex] > par_vertex_of_simpl)
                par_vertex_of_simpl = map_par_vertices[vertex];
        }

        vector<double>::iterator lower1;
        lower1 = lower_bound(bins1.begin(), bins1.end(), par_vertex_of_simpl);
        size_t index_vertex_in_bins1(lower1 - bins1.begin());
        map_simplices[index_vertex_in_bins1].push_back(index_simpl);
    }

    // loop on indices in bins1
    // set the cumulative sum of changes equal to rows in surface
    vector<int> changes(nbins2+1, 0);
    vector<vector<int>> surface(nbins1, vector< int >(nbins2, 0));
    vector<int> possible_changes{1, -1, 1};

    for (size_t ind_bins1(0); ind_bins1 < bins1.size(); ind_bins1++)
    {
        vector<size_t> args_simplices_to_insert(map_simplices[ind_bins1]);

        // loop to update changes
        for (const size_t &arg: args_simplices_to_insert)
        {
            vector<int> inserted_simpl(sorted_simplices[arg]);
            double inserted_par(sorted_par_simplices[arg]);

            size_t dim(dim_simplex(inserted_simpl));

            vector<double>::iterator lower2;
            lower2 = lower_bound(bins2.begin(), bins2.end(), inserted_par);
            changes[(lower2 - bins2.begin())] += possible_changes[dim];
        }
        // cumulative sum into surface
        if (ind_bins1 != 0)
        {
            int c(changes[0]);
            for (size_t index(0); index < nbins2; index++)
            {
                surface[ind_bins1-1][index] = c + changes[index+1];
                c += changes[index+1];
            }
        }

    }

    return surface;
}

//================================================

vector<vector<int>> bifiltration_3d(vector<double> unsorted_par_vertices,
                                    vector<double> sorted_par_simplices,
                                    vector<vector<int>> sorted_simplices,
                                    size_t nbins1, size_t nbins2,
                                    vector<double> minmax1,
                                    vector<double> minmax2)
{
    // sorted indices of vertices
    vector<size_t> argvert(unsorted_par_vertices.size());
    iota(argvert.begin(), argvert.end(), 0);

    // map of vertices parametrization
    map< size_t, double > map_par_vertices;
    for (size_t k(0); k < unsorted_par_vertices.size(); k++)
        map_par_vertices[k] = unsorted_par_vertices[k];

    // sort indexes based on comparing values in unsorted_par_vertices
    sort(argvert.begin(), argvert.end(),
        [&unsorted_par_vertices](size_t i1, size_t i2)
            {return unsorted_par_vertices[i1] < unsorted_par_vertices[i2];});

    // make bins1-->(vertices) and bins2-->(simplices)

    double step1((minmax1[1] - minmax1[0]) / nbins1);
    double step2((minmax2[1] - minmax2[0]) / nbins2);

    vector<double> bins1(nbins1+1);   // the last element of these
    vector<double> bins2(nbins2+1);   // is different from np.linspace()

    for (size_t row(0); row < nbins1+1; row++)
        bins1[row] = step1 * row + minmax1[0];

    for (size_t col(0); col < nbins2+1; col++)
        bins2[col] = step2 * col + minmax2[0];

    // map of simplices representing bi-filtration
    map<size_t, vector<size_t>> map_simplices;
    for (size_t index_simpl(0); index_simpl < sorted_simplices.size(); index_simpl++)
    {
        vector<int> simpl(sorted_simplices[index_simpl]);

        double par_vertex_of_simpl(map_par_vertices[simpl[0]]);
        for (const int &vertex: simpl)
        {
            if (vertex == -1)
                continue;  // if simp is "sorted" use break

            if (map_par_vertices[vertex] > par_vertex_of_simpl)
                par_vertex_of_simpl = map_par_vertices[vertex];
        }

        vector<double>::iterator lower1;
        lower1 = lower_bound(bins1.begin(), bins1.end(), par_vertex_of_simpl);
        size_t index_vertex_in_bins1(lower1 - bins1.begin());
        map_simplices[index_vertex_in_bins1].push_back(index_simpl);
    }

    // loop on indices in bins1
    // set the cumulative sum of changes equal to rows in surface
    vector<int> changes(nbins2+1, 0);
    vector<vector<int>> surface(nbins1, vector< int >(nbins2, 0));
    vector<int> possible_changes{1, -1, 1, -1};

    for (size_t ind_bins1(0); ind_bins1 < bins1.size(); ind_bins1++)
    {
        vector<size_t> args_simplices_to_insert(map_simplices[ind_bins1]);

        // loop to update changes
        for (const size_t &arg: args_simplices_to_insert)
        {
            vector<int> inserted_simpl(sorted_simplices[arg]);
            double inserted_par(sorted_par_simplices[arg]);

            size_t dim(dim_simplex(inserted_simpl));

            vector<double>::iterator lower2;
            lower2 = lower_bound(bins2.begin(), bins2.end(), inserted_par);
            changes[(lower2 - bins2.begin())] += possible_changes[dim];
        }
        // cumulative sum into surface
        if (ind_bins1 != 0)
        {
            int c(changes[0]);
            for (size_t index(0); index < nbins2; index++)
            {
                surface[ind_bins1-1][index] = c + changes[index+1];
                c += changes[index+1];
            }
        }

    }

    return surface;
}

