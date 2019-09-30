#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>       // std::vector
#include <numeric>      // std::partial_sum
#include "libutils.hpp"

namespace py = pybind11;

//================================================

std::vector<std::vector<int>> pad(const std::vector<std::vector<int>> &image)
{
    size_t num_rows = image.size();
    size_t num_cols = image[0].size();

    std::vector<std::vector<int>> padded(num_rows+2, std::vector<int>(num_cols+2, 255));

    for (size_t i = 0; i < num_rows; i++)
    {
        for (size_t j = 0; j < num_cols; j++)
        {
            padded[i+1][j+1] = image[i][j];
        }
    }

    return padded;
}

//================================================

size_t num_after(const std::vector<std::vector<int>> &padded, const std::vector<std::vector<size_t>> &indices)
{
    size_t num = 0;
    int pixel_value = padded[indices[4][0]][indices[4][1]];
    const std::vector<size_t> powers_of_two = {1, 2, 4, 8, 0,
                                                16, 32, 64, 128};    
    // Before and including central_pixel
    for (size_t k = 0; k < 5; k++)
    {
        if ( padded[indices[k][0]][indices[k][1]] <= pixel_value )
            num |= powers_of_two[k];
    }
    // After central_pixel
    for (size_t k = 5; k < 9; k++)
    {
        if ( padded[indices[k][0]][indices[k][1]] < pixel_value )
            num |= powers_of_two[k];
    }

    return num;
}

//================================================

std::vector<int> naive_image_2d(std::vector<std::vector<int>> image, int M)
{
    std::vector< int > ecc(M, 0);

    for (size_t i = 0; i < M; i++)
    {
      std::vector< std::vector< bool > > thresh = threshold_image_2d(image, static_cast<int>(i));
        ecc[i] = char_binary_image_2d(thresh);
    }

    return ecc;
}

//================================================

std::vector<int> image_2d(const std::vector<std::vector<int>> & image,
                          int M,
			  const std::vector< int > &vec_char)
{
    size_t num_rows = image.size();
    size_t num_cols = image[0].size();

    // Euler characteristic curve
    std::vector< int > ecc(M, 0);

    // Padded image
    std::vector<std::vector<int>> padded(num_rows+2, std::vector<int>(num_cols+2, M-1));
    for (size_t i = 0; i < num_rows; i++)
    {
        for (size_t j = 0; j < num_cols; j++)
        {
            padded[i+1][j+1] = image[i][j];
        }
    }

    // Loop over all pixel in original image
    for (size_t i = 1; i < num_rows+1; i++)
    {
        for (size_t j = 1; j < num_cols+1; j++)
        {
            int pixel_value = padded[i][j];

            std::vector<std::vector<size_t>> indices(9, std::vector<size_t>(2));
            indices[0] = {i-1, j-1};
            indices[1] = {i-1, j};
            indices[2] = {i-1, j+1};
            indices[3] = {i, j-1};
            indices[4] = {i, j};
            indices[5] = {i, j+1};
            indices[6] = {i+1, j-1};
            indices[7] = {i+1, j};
            indices[8] = {i+1, j+1};

            size_t num = num_after(padded, indices);

            // Add euler change at pixel_value
            ecc[static_cast<size_t>(pixel_value)] += static_cast<int>(vec_char[num]);
        }
    }

    // Cumulative sum of euler changes
    std::partial_sum(ecc.begin(), ecc.end(), ecc.begin());

    return ecc;
}
