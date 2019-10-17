#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <numeric>              // std::partial_sum, std::iota
#include <algorithm>            // std::copy, std::lower_bound

using namespace std;

//================================================

int sum_bool_2d(vector<vector<bool>> matrix)
{
    size_t numI = matrix.size();
    size_t numJ = matrix[0].size();
    
    int total = 0;
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            if (matrix[i][j] == true)
                total += 1;
        }
    }

    return total;
}

//================================================

vector<vector<int>> pad_2d(const vector<vector<int>> &image, int M)
{
    size_t numI = image.size();
    size_t numJ = image[0].size();
    
    vector<vector<int>> padded(numI+2, vector<int>(numJ+2, M));

    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            padded[i+1][j+1] = image[i][j];
        }
    }

    return padded;
}

//================================================

vector<vector<bool>> threshold_image_2d(const vector<vector<int>> &image, int value)
{
    size_t numI = image.size();
    size_t numJ = image[0].size();

    vector<vector<bool>> binary_thresh(numI, vector<bool>(numJ, false));

    // Loop over all pixel in original image
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            if (image[i][j] <= value)
            {
                binary_thresh[i][j] = true;
            }
        }
    }

    return binary_thresh;
}

//================================================

vector<vector<bool>> elementwise_AND_2d(const vector<vector<bool>> &image1, const vector<vector<bool>> &image2)
{
    size_t numI = image1.size();
    size_t numJ = image1[0].size();

    vector<vector<bool>> result(numI, vector<bool>(numJ, false));

    // Loop over all pixel in original image
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
	  if (image1[i][j] == true && image2[i][j] == true) {
                result[i][j] = true;
            }
        }
    }

    return result;
}

//================================================

vector<vector<int>> neigh_pixel_2d(const vector<vector<int>> &padded, size_t i , size_t j)
{
    vector<vector<int>> neigh_3_3(3, vector<int>(3, 0));

    neigh_3_3[0][0] = padded[i-1][j-1];
    neigh_3_3[1][0] = padded[i][j-1];
    neigh_3_3[2][0] = padded[i+1][j-1];
    neigh_3_3[0][1] = padded[i-1][j];
    neigh_3_3[1][1] = padded[i][j];
    neigh_3_3[2][1] = padded[i+1][j];
    neigh_3_3[0][2] = padded[i-1][j+1];
    neigh_3_3[1][2] = padded[i][j+1];
    neigh_3_3[2][2] = padded[i+1][j+1];

    return neigh_3_3;
}

//================================================

vector<vector<bool>> binary_neigh_pixel_2d(const vector<vector<int>> &padded, size_t i , size_t j , int pixel_value)
{
    vector<vector<int>> neigh_3_3 = neigh_pixel_2d(padded, i, j);
    
    vector<vector<bool>> binary_neigh_3_3(3, vector<bool>(3, false));

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            // need to take account of order pixels
            // strinct inequality after central pixel
            size_t k = j + i*3;
            
            if (k < 5) {  // up to central pixel included
                if (neigh_3_3[i][j] <= pixel_value)
                    binary_neigh_3_3[i][j] = true;
            } else {
                if (neigh_3_3[i][j] < pixel_value)
                    binary_neigh_3_3[i][j] = true;
            }
        }
    }

    return binary_neigh_3_3;
}

//================================================

int char_binary_image_2d(vector<vector<bool>> input)
{
    // Input shape: binary image number of rows and columns
    size_t numI = input.size();
    size_t numJ = input[0].size();

    // Matrices for vectices, horizontal edges and vertical edges
    vector<vector<bool>> V(numI+1, vector<bool>(numJ+1, false));
    vector<vector<bool>> Eh(numI+1, vector<bool>(numJ, false));
    vector<vector<bool>> Ev(numI, vector<bool>(numJ+1, false));

    // Loop over pixels to update V, Eh, Ev
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            if (input[i][j] == true) {
                V[i][j]     = true;
                V[i+1][j]   = true;
                V[i][j+1]   = true;
                V[i+1][j+1] = true;

                Eh[i][j]    = true;
                Eh[i+1][j]  = true;

                Ev[i][j]    = true;
                Ev[i][j+1]  = true;
            }
        }
    }

    // Sum of elements in the matrices
    int v  = sum_bool_2d(V);
    int eh = sum_bool_2d(Eh);
    int ev = sum_bool_2d(Ev);
    int f  = sum_bool_2d(input);

    int EC = v - eh - ev + f;
    return EC;
}

//================================================

size_t number_from_neigh_2d(vector<vector<bool>> neigh)
{
    const vector<size_t> powers_of_two = {1, 2, 4, 8, 0,
                                          16, 32, 64, 128};
    size_t num = 0;

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            if (neigh[i][j] == true) {
                num |= powers_of_two[j+i*3];
            }
        }
    }

    return num;
}

//================================================

int sum_bool_3d(vector<vector<vector<bool>>> matrix)
{
    size_t numI = matrix.size();
    size_t numJ = matrix[0].size();
    size_t numK = matrix[0][0].size();

    int total = 0;
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            for (size_t k = 0; k < numK; ++k) {
                if (matrix[i][j][k] == true)
                    total += 1;
            }
        }
    }
    
    return total;
}

//================================================

vector<vector<vector<int>>> pad_3d(const vector<vector<vector<int>>> &image, int M)
{
    size_t numI = image.size();
    size_t numJ = image[0].size();
    size_t numK = image[0][0].size();
    
    vector<vector<vector<int>>> padded(numI+2, vector<vector<int>>(numJ+2, vector<int>(numK+2, M)));

    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            for (size_t k = 0; k < numK; ++k) {
                padded[i+1][j+1][k+1] = image[i][j][k];
            }
        }
    }

    return padded;
}

//================================================

vector<vector<vector<bool>>> threshold_image_3d(const vector<vector<vector<int>>> &image, int value) 
{
    size_t numI = image.size();
    size_t numJ = image[0].size();
    size_t numK = image[0][0].size();

    vector<vector<vector<bool>>> binary_thresh(numI, vector<vector<bool>>(numJ, vector<bool>(numK, false)));

    // Loop over all pixel in original image
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            for (size_t k = 0; k < numK; ++k) {
                if (image[i][j][k] <= value) {
                    binary_thresh[i][j][k] = true;
                }
            }
        }
    }

    return binary_thresh;
}

//================================================

vector<vector<vector<bool>>> elementwise_AND_3d(const vector<vector<vector<bool>>> &image1, const vector<vector<vector<bool>>> &image2)
{
    size_t numI = image1.size();
    size_t numJ = image1[0].size();
    size_t numK = image1[0][0].size();

    vector<vector<vector<bool>>> result(numI, vector<vector<bool>>(numJ, vector<bool>(numK, false)));

    // Loop over all pixel in original image
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            for (size_t k = 0; k < numK; ++k) {
                if (image1[i][j][k] == true && image2[i][j][k] == true) {
                    result[i][j][k] = true;
                }

            }
        }

    }
    
    return result;
}

//================================================

vector<vector<vector<int>>> neigh_voxel_3d(const vector<vector<vector<int>>> &padded, size_t i , size_t j , size_t k)
{
    vector<vector<vector<int>>> neigh_3_3_3(3, vector<vector<int>>(3, vector<int>(3, 0)));

    neigh_3_3_3[0][0][0] = padded[i-1][j-1][k-1];
    neigh_3_3_3[1][0][0] = padded[i][j-1][k-1];
    neigh_3_3_3[2][0][0] = padded[i+1][j-1][k-1];
    neigh_3_3_3[0][1][0] = padded[i-1][j][k-1];
    neigh_3_3_3[1][1][0] = padded[i][j][k-1];
    neigh_3_3_3[2][1][0] = padded[i+1][j][k-1];
    neigh_3_3_3[0][2][0] = padded[i-1][j+1][k-1];
    neigh_3_3_3[1][2][0] = padded[i][j+1][k-1];
    neigh_3_3_3[2][2][0] = padded[i+1][j+1][k-1];

    neigh_3_3_3[0][0][1] = padded[i-1][j-1][k];
    neigh_3_3_3[1][0][1] = padded[i][j-1][k];
    neigh_3_3_3[2][0][1] = padded[i+1][j-1][k];
    neigh_3_3_3[0][1][1] = padded[i-1][j][k];
    neigh_3_3_3[1][1][1] = padded[i][j][k];
    neigh_3_3_3[2][1][1] = padded[i+1][j][k];
    neigh_3_3_3[0][2][1] = padded[i-1][j+1][k];
    neigh_3_3_3[1][2][1] = padded[i][j+1][k];
    neigh_3_3_3[2][2][1] = padded[i+1][j+1][k];

    neigh_3_3_3[0][0][2] = padded[i-1][j-1][k+1];
    neigh_3_3_3[1][0][2] = padded[i][j-1][k+1];
    neigh_3_3_3[2][0][2] = padded[i+1][j-1][k+1];
    neigh_3_3_3[0][1][2] = padded[i-1][j][k+1];
    neigh_3_3_3[1][1][2] = padded[i][j][k+1];
    neigh_3_3_3[2][1][2] = padded[i+1][j][k+1];
    neigh_3_3_3[0][2][2] = padded[i-1][j+1][k+1];
    neigh_3_3_3[1][2][2] = padded[i][j+1][k+1];
    neigh_3_3_3[2][2][2] = padded[i+1][j+1][k+1];

    return neigh_3_3_3;
}

//================================================

vector<vector<vector<bool>>> binary_neigh_voxel_3d(const vector<vector<vector<int>>> &padded, size_t i , size_t j , size_t k, int voxel_value)
{
    vector<vector<vector<int>>> neigh_3_3_3 = neigh_voxel_3d(padded, i, j, k);

    vector<vector<vector<bool>>> binary_neigh_3_3_3(3, vector<vector<bool>>(3, vector<bool>(3, false)));

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                // need to take account of order pixels
                // strinct inequality after central pixel
                size_t s = k + j*3 + i*9;
            
                if (s < 14) {  // up to central pixel included
                    if (neigh_3_3_3[i][j][k] <= voxel_value)
                        binary_neigh_3_3_3[i][j][k] = true;
                } else {
                    if (neigh_3_3_3[i][j][k] < voxel_value)
                        binary_neigh_3_3_3[i][j][k] = true;
                }   
            }
        }
    }
    

    return binary_neigh_3_3_3;
}

//================================================

int char_binary_image_3d(vector<vector<vector<bool>>> input)
{
    size_t numI = input.size();
    size_t numJ = input[0].size();
    size_t numK = input[0][0].size();

    vector<vector<vector<bool>>> V(numI+1, vector<vector<bool>>(numJ+1, vector<bool>(numK+1, false)));
    vector<vector<vector<bool>>> Ei(numI, vector<vector<bool>>(numJ+1, vector<bool>(numK+1, false)));
    vector<vector<vector<bool>>> Ej(numI+1, vector<vector<bool>>(numJ, vector<bool>(numK+1, false)));
    vector<vector<vector<bool>>> Ek(numI+1, vector<vector<bool>>(numJ+1, vector<bool>(numK, false)));
    vector<vector<vector<bool>>> Fij(numI, vector<vector<bool>>(numJ, vector<bool>(numK+1, false)));
    vector<vector<vector<bool>>> Fik(numI, vector<vector<bool>>(numJ+1, vector<bool>(numK, false)));
    vector<vector<vector<bool>>> Fjk(numI+1, vector<vector<bool>>(numJ, vector<bool>(numK, false)));
    vector<vector<vector<bool>>> C(numI, vector<vector<bool>>(numJ, vector<bool>(numK, false)));
        
    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            for (size_t k = 0; k < numK; ++k) {
                if (input[i][j][k] == true) {
                    V[i][j][k]       = true;
                    V[i+1][j][k]     = true;
                    V[i][j+1][k]     = true;
                    V[i][j][k+1]     = true;
                    V[i+1][j+1][k]   = true;
                    V[i+1][j][k+1]   = true;
                    V[i][j+1][k+1]   = true;
                    V[i+1][j+1][k+1] = true;
                    
                    Ei[i][j][k]      = true;
                    Ei[i][j+1][k]    = true;
                    Ei[i][j][k+1]    = true;
                    Ei[i][j+1][k+1]  = true;
                    
                    Ej[i][j][k]      = true;
                    Ej[i+1][j][k]    = true;
                    Ej[i][j][k+1]    = true;
                    Ej[i+1][j][k+1]  = true;
                    
                    Ek[i][j][k]      = true;
                    Ek[i][j+1][k]    = true;
                    Ek[i+1][j][k]    = true;
                    Ek[i+1][j+1][k]  = true;
                    
                    Fij[i][j][k]     = true;
                    Fij[i][j][k+1]   = true;
                    
                    Fik[i][j][k]     = true;
                    Fik[i][j+1][k]   = true;
                    
                    Fjk[i][j][k]     = true;
                    Fjk[i+1][j][k]   = true;
                    
                    C[i][j][k]       = true;
                }
            }
        }
    }

    int v  = sum_bool_3d(V);
    int ei = sum_bool_3d(Ei);
    int ej = sum_bool_3d(Ej);
    int ek = sum_bool_3d(Ek);
    int fij = sum_bool_3d(Fij);
    int fik = sum_bool_3d(Fik);
    int fjk = sum_bool_3d(Fjk);
    int c = sum_bool_3d(C);

    int EC = v - ei - ej - ek + fij + fik + fjk - c;
    
    return EC;
}

//================================================

size_t number_from_neigh_3d(vector<vector<vector<bool>>> neigh)
{
    size_t num = 0;

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j){
            for (size_t k = 0; k < 3; ++k) {
                size_t bits_shift = static_cast<size_t>(1ULL << (k + j*3 + i * 9));
                size_t bits_shift_minus_1 = static_cast<size_t>(1ULL << (k + j*3 + i * 9 - 1));
                if (bits_shift < 8192) { // 2^13 == 8192, and 14th
                                         // voxel is central one
                    if (neigh[i][j][k] == true)
                        num |= bits_shift;
                } else if (bits_shift > 8192)
                {
                    if (neigh[i][j][k] == true)
                        num |= bits_shift_minus_1;
                }
            }
        }
    }

    return num;
}

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

vector<int> image_3d(const vector<vector<vector<int>>> &image,
                     const vector<int> &vector_euler_changes,
                     int M)
{
    size_t numI = image.size();
    size_t numJ = image[0].size();
    size_t numK = image[0][0].size();

    vector<int> ecc(M+1, 0);

    vector<vector<vector<int>>> padded = pad_3d(image, M);

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

vector<int> filtration(vector<int>    dim_simplices,
                       vector<double> parametrization,
                       vector<double> bins)
{
    size_t num_elements = bins.size();
   
    vector<int> euler_char_curve(num_elements, 0);
  
    // possible changes due to addition of vertex edge or triangle
    vector<int> possible_changes{1, -1, 1, -1};

    // loop on simplices and update euler curve
    for (size_t i = 0; i < dim_simplices.size(); ++i) {
        size_t dim_simplex = dim_simplices[i];
        double par = parametrization[i];
        
        vector<double>::iterator lower;
        lower = lower_bound(bins.begin(), bins.end(), par);
        euler_char_curve[(lower - bins.begin())] += possible_changes[dim_simplex];
    }

    int tmp = euler_char_curve[0];
    for (size_t s = 1; s < euler_char_curve.size(); ++s) {
        euler_char_curve[s] += tmp;
        tmp = euler_char_curve[s];
    }
    
    return euler_char_curve;
}

//================================================

PYBIND11_MODULE(curve, m) {
    m.doc() = "Euler characteristic curves cpp bindings.";

    m.def("naive_image_2d", &naive_image_2d);
    m.def("image_2d", &image_2d);
    m.def("naive_image_3d", &naive_image_3d);
    m.def("image_3d", &image_3d);
    m.def("filtration", &filtration);
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
