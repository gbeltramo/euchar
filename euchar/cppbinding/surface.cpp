#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>            // std::copy, std::lower_bound
#include <vector>

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

vector<vector<int>> naive_images_2d(const vector<vector<int>> &image1,
                                    const vector<vector<int>> &image2,
                                    int M1, int M2)
{
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

vector<vector<int>> images_3d(const vector<vector<vector<int>>> &image1,
			      const vector<vector<vector<int>>> &image2,
                              vector<int> vector_euler_changes,
                              int M1 = 255, int M2 = 255)
{
    size_t numI = image1.size();
    size_t numJ = image1[0].size();
    size_t numK = image1[0][0].size();
    
    vector<vector<int>> surface(M1+1, vector<int>(M2+1, 0));
    
    // Padded images
    vector<vector<vector<int>>> padded1 = pad_3d(image1, M1);
    vector<vector<vector<int>>> padded2 = pad_3d(image2, M2);

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

//================================================

PYBIND11_MODULE(surface, m) {
    m.doc() = "Euler characteristic surfaces cpp bindings.";

    m.def("naive_images_2d", &naive_images_2d);
    m.def("images_2d", &images_2d);
    m.def("naive_images_3d", &naive_images_3d);
    m.def("images_3d", &images_3d);
    m.def("bifiltration", &bifiltration);
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
