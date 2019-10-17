#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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

PYBIND11_MODULE(utils, m) {
    m.doc() = "";

    // 2D functions
    m.def("sum_bool_2d", &sum_bool_2d);
    m.def("char_binary_image_2d", &char_binary_image_2d);
    m.def("neigh_2d_from_number", &neigh_2d_from_number);
    m.def("number_from_neigh_2d", &number_from_neigh_2d);
    m.def("vector_of_euler_changes_2d", &vector_of_euler_changes_2d);

    // 3D functions
    m.def("sum_bool_3d", &sum_bool_3d);
    m.def("char_binary_image_3d", &char_binary_image_3d);
    m.def("neigh_3d_from_number", &neigh_3d_from_number);
    m.def("number_from_neigh_3d", &number_from_neigh_3d);
    m.def("vector_of_euler_changes_3d", &vector_of_euler_changes_3d);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
