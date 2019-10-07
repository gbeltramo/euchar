#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;
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

vector<vector<bool>> neigh_2d_from_number(size_t num)
{
  const vector<size_t> powers_of_two = {1, 2, 4, 8, 0,
                                              16, 32, 64, 128};
  vector<vector<bool>> neigh(3, vector<bool>(3, false));
  for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            if (num & powers_of_two[j+i*3] || (i == 1 && j == 1))
                neigh[i][j] = true;
        }
    }

    return neigh;
}

//================================================

size_t number_from_neigh_2d(vector<vector<bool>> neigh)
{
    const vector<size_t> powers_of_two = {1, 2, 4, 8, 0,
                                          16, 32, 64, 128};

    if ( neigh.size() != 3 || neigh[0].size() != 3)
        throw runtime_error("Nieghbor must be 3x3 matrix");

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

vector<short> vector_of_euler_changes_2d()
{
    vector<short> euler(256, 0);

    for (size_t num = 0; num < 256; ++num) {
        vector<vector<bool>> binary_neigh = neigh_2d_from_number(num);
        int ec_after = char_binary_image_2d(binary_neigh);
        binary_neigh[1][1] = false;
        int ec_before = char_binary_image_2d(binary_neigh);
        euler[num] = static_cast<short>(ec_after - ec_before);
      }

    return euler;
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

vector<vector<vector<int>>> pad_3d(py::array_t<int> input, int M)
{
    auto image = input.unchecked<3>();
    
    size_t numI = image.shape(0);
    size_t numJ = image.shape(1);
    size_t numK = image.shape(2);
    
    vector<vector<vector<int>>> padded(numI+2, vector<vector<int>>(numJ+2, vector<int>(numK+2, M)));

    for (size_t i = 0; i < numI; ++i) {
        for (size_t j = 0; j < numJ; ++j) {
            for (size_t k = 0; k < numK; ++k) {
                padded[i+1][j+1][k+1] = image(i, j, k);
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

vector<vector<vector<bool>>> neigh_3d_from_number(size_t num)
{
    vector<vector<vector<bool>>> neigh(3, vector<vector<bool>>(3, vector<bool>(3, false)));

    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                size_t bits_shift = static_cast<size_t>(1ULL << (k + j*3 + i * 9));
                size_t bits_shift_minus_1 = static_cast<size_t>(1ULL << (k + j*3 + i * 9 - 1));
                if (bits_shift < 8192)
                {
                    if (num & bits_shift)
                        neigh[i][j][k] = true;
                } else if (bits_shift == 8192)
                {
                    neigh[i][j][k] = true;
                } else
                {
                    if (num & bits_shift_minus_1)
                        neigh[i][j][k] = true;
                }
                
            }
        }
    }

    return neigh;
    
}

//================================================

size_t number_from_neigh_3d(vector<vector<vector<bool>>> neigh)
{
    if ( neigh.size() != 3 || neigh[0].size() != 3 || neigh[0][0].size() != 3 )
        throw runtime_error("Nieghbor must be 3x3x3 matrix");

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

vector<short> vector_of_euler_changes_3d(size_t offset=0, size_t max_value=67108864)
{
    vector<short> euler(max_value, 0);

    for (size_t num = offset; num < euler.size(); num++)
    {
        vector<vector<vector<bool>>> binary_neigh = neigh_3d_from_number(static_cast<size_t>(num));
        
        int ec_after = char_binary_image_3d(binary_neigh);
        binary_neigh[1][1][1] = false;
        int ec_before = char_binary_image_3d(binary_neigh);
        euler[num] = static_cast<short>(ec_after - ec_before);
      }

  return euler;
}

//================================================

size_t dim_simplex(const vector<int> &simplex)
{
    size_t length = simplex.size();
    size_t equal_minus_one(0);
    
    for (const int &ind: simplex) {
        if (ind == -1)
            equal_minus_one++;
    }

    return length - equal_minus_one - 1;
}

//================================================

vector<bool> filter_parametrization_on(const vector<int> &vertices,
                                       const vector<vector<int>> &simplices,
                                       const vector<double> &param)
{
    size_t L = simplices.size();
    vector<bool> output(L, true);

    // make vertices into set for fast lookup
    set<int> set_vertices;
    copy(vertices.begin(),
              vertices.end(),
              inserter(set_vertices, set_vertices.end()));

    // Loop on elements of simplices and update mask
    for (size_t k = 0; k < L; k++) {
        vector<int> simpk(simplices[k]);
        size_t dim(dim_simplex(simpk));

        bool simplex_in(true);
        for (size_t i=0; i < dim+1; i++) {
            if (set_vertices.find(simpk[i]) == set_vertices.end()) {
                simplex_in = false;
                break;
            }
        }

        if (simplex_in == false) {
            output[k] = false;
        }
    }

    return output;
}

