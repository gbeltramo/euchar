#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

//================================================

int sum_bool_matrix(std::vector<std::vector<bool>> matrix)
{
    size_t num_rows = matrix.size();
    size_t num_cols = matrix[0].size();

    int total = 0;
    for (size_t i = 0; i < num_rows; i++)
    {
        for (size_t j = 0; j < num_cols; j++)
        {
            if (matrix[i][j] == true)
                total += 1;
        }
    }

    return total;
}

//================================================

std::vector<std::vector<bool>> threshold_image_2d(const std::vector<std::vector<int>> &image, int value)
{
    size_t num_rows = image.size();
    size_t num_cols = image[0].size();

    std::vector<std::vector<bool>> binary_thresh(num_rows, std::vector<bool>(num_cols, false));

    // Loop over all pixel in original image
    for (size_t i = 0; i < num_rows; i++)
    {
        for (size_t j = 0; j < num_cols; j++)
        {
            if (image[i][j] <= value)
            {
                binary_thresh[i][j] = true;
            }
        }
    }

    return binary_thresh;
}

//================================================

std::vector<std::vector<bool>> and_two_binary_images(const std::vector<std::vector<bool>> &binary1, const std::vector<std::vector<bool>> &binary2)
{
    size_t num_rows = binary1.size();
    size_t num_cols = binary1[0].size();

    std::vector<std::vector<bool>> result(num_rows, std::vector<bool>(num_cols, false));

    // Loop over all pixel in original image
    for (size_t i = 0; i < num_rows; i++)
    {
        for (size_t j = 0; j < num_cols; j++)
        {
	  if (binary1[i][j] == true && binary2[i][j] == true)
            {
                result[i][j] = true;
            }
        }
    }

    return result;
}

//================================================

int char_binary_image_2d(std::vector<std::vector<bool>> input)
{
    // Input shape: binary image number of rows and columns
    size_t num_rows = input.size();
    size_t num_cols = input[0].size();

    // Matrices for vectices, horizontal edges and vertical edges
    std::vector<std::vector<bool>> V(num_rows+1, std::vector<bool>(num_cols+1, false));
    std::vector<std::vector<bool>> Eh(num_rows+1, std::vector<bool>(num_cols, false));
    std::vector<std::vector<bool>> Ev(num_rows, std::vector<bool>(num_cols+1, false));

    // Loop over pixels to update V, Eh, Ev
    for (size_t i = 0; i < num_rows; i++)
    {
        for (size_t j = 0; j < num_cols; j++)
        {
            if (input[i][j] == true)
            {
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
    int v  = sum_bool_matrix(V);
    int eh = sum_bool_matrix(Eh);
    int ev = sum_bool_matrix(Ev);
    int f  = sum_bool_matrix(input);

    int EC = v - eh - ev + f;
    return EC;
}

//================================================

std::vector<std::vector<bool>> neigh_2d_from_number(uint8_t num)
{
  const std::vector<uint8_t> powers_of_two = {1, 2, 4, 8, 0,
                                              16, 32, 64, 128};
  std::vector<std::vector<bool>> neigh(3, std::vector<bool>(3, false));
  for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            if (num & powers_of_two[j+i*3] || (i == 1 && j == 1))
                neigh[i][j] = true;
        }

    }

    return neigh;
}

//================================================

uint8_t number_from_neigh_2d(std::vector<std::vector<bool>> neigh)
{
    const std::vector<uint8_t> powers_of_two = {1, 2, 4, 8, 0,
                                              16, 32, 64, 128};

    if ( neigh.size() != 3 || neigh[0].size() != 3)
        throw std::runtime_error("Nieghbor must be 3x3 matrix");

    uint8_t num = 0;

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            if (neigh[i][j] == true)
            {
                num |= powers_of_two[j+i*3];
            }
        }
    }

    return num;
}

//================================================

std::vector<short> vector_of_euler_changes_2d()
{
    std::vector<short> euler(256, 0);

    for (size_t num = 0; num < 256; num++)
      {
          std::vector<std::vector<bool>> binary_neigh = neigh_2d_from_number(static_cast<uint8_t>(num));
        int ec_after = char_binary_image_2d(binary_neigh);
        binary_neigh[1][1] = false;
        int ec_before = char_binary_image_2d(binary_neigh);
        euler[num] = static_cast<short>(ec_after - ec_before);
      }

    return euler;
}

//================================================

int sum_bool_matrix_3d(std::vector<std::vector<std::vector<bool>>> matrix)
{
    size_t numI = matrix.size();
    size_t numJ = matrix[0].size();
    size_t numK = matrix[0][0].size();

    int total = 0;
    for (size_t i = 0; i < numI; i++)
    {
        for (size_t j = 0; j < numJ; j++)
        {
            for (size_t k = 0; k < numK; k++)
            {
                if (matrix[i][j][k] == true)
                    total += 1;
            }
        }
    }
    
    return total;
}

//================================================

int char_binary_image_3d(std::vector<std::vector<std::vector<bool>>> input)
{
    size_t numI = input.size();
    size_t numJ = input[0].size();
    size_t numK = input[0][0].size();

    std::vector<std::vector<std::vector<bool>>> V(numI+1, std::vector<std::vector<bool>>(numJ+1, std::vector<bool>(numK+1, false)));
    std::vector<std::vector<std::vector<bool>>> Ei(numI, std::vector<std::vector<bool>>(numJ+1, std::vector<bool>(numK+1, false)));
    std::vector<std::vector<std::vector<bool>>> Ej(numI+1, std::vector<std::vector<bool>>(numJ, std::vector<bool>(numK+1, false)));
    std::vector<std::vector<std::vector<bool>>> Ek(numI+1, std::vector<std::vector<bool>>(numJ+1, std::vector<bool>(numK, false)));
    std::vector<std::vector<std::vector<bool>>> Fij(numI, std::vector<std::vector<bool>>(numJ, std::vector<bool>(numK+1, false)));
    std::vector<std::vector<std::vector<bool>>> Fik(numI, std::vector<std::vector<bool>>(numJ+1, std::vector<bool>(numK, false)));
    std::vector<std::vector<std::vector<bool>>> Fjk(numI+1, std::vector<std::vector<bool>>(numJ, std::vector<bool>(numK, false)));
    std::vector<std::vector<std::vector<bool>>> C(numI, std::vector<std::vector<bool>>(numJ, std::vector<bool>(numK, false)));
        
    for (size_t i = 0; i < numI; i++)
    {
        for (size_t j = 0; j < numJ; j++)
        {
            for (size_t k = 0; k < numK; k++)
            {
                if (input[i][j][k] == true)
                {
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

    int v  = sum_bool_matrix_3d(V);
    int ei = sum_bool_matrix_3d(Ei);
    int ej = sum_bool_matrix_3d(Ej);
    int ek = sum_bool_matrix_3d(Ek);
    int fij = sum_bool_matrix_3d(Fij);
    int fik = sum_bool_matrix_3d(Fik);
    int fjk = sum_bool_matrix_3d(Fjk);
    int c = sum_bool_matrix_3d(C);

    int EC = v - ei - ej - ek + fij + fik + fjk - c;
    
    return EC;
}

//================================================

std::vector<std::vector<std::vector<bool>>> neigh_3d_from_number(size_t num)
{
    std::vector<std::vector<std::vector<bool>>> neigh(3, std::vector<std::vector<bool>>(3, std::vector<bool>(3, false)));

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            for (size_t k = 0; k < 3; k++)
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

size_t number_from_neigh_3d(std::vector<std::vector<std::vector<bool>>> neigh)
{
    if ( neigh.size() != 3 || neigh[0].size() != 3 || neigh[1].size() != 3 )
        throw std::runtime_error("Nieghbor must be 3x3x3 matrix");

    size_t num = 0;

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                size_t bits_shift = static_cast<size_t>(1ULL << (k + j*3 + i * 9));
                size_t bits_shift_minus_1 = static_cast<size_t>(1ULL << (k + j*3 + i * 9 - 1));
                if (bits_shift < 8192)
                {
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

std::vector<short> vector_of_euler_changes_3d(size_t offset=0, size_t max_value=6710886)
{
    std::vector<short> euler(max_value, 0);

    for (size_t num = offset; num < euler.size(); num++)
    {
        std::vector<std::vector<std::vector<bool>>> binary_neigh = neigh_3d_from_number(static_cast<size_t>(num));
        
        int ec_after = char_binary_image_3d(binary_neigh);
        binary_neigh[1][1][1] = false;
        int ec_before = char_binary_image_3d(binary_neigh);
        euler[num] = static_cast<short>(ec_after - ec_before);
      }

  return euler;
}

  
