#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "libutils.hpp"
#include "libcurve.hpp"

//================================================

std::vector<std::vector<int>> naive_surface_2d(const std::vector<std::vector<int>> &image1,
    const std::vector<std::vector<int>> &image2,
    int M)
{
    // Check length of key
    if ( image1.size() != image2.size() || image1[0].size() != image2[0].size())
        throw std::runtime_error("Images dimensions must be equal.");

    // Output surface
    std::vector<std::vector<int>> surface(M, std::vector<int>(M, 0));

    // Loop over all pixel in original image
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            std::vector<std::vector<bool>> thresh1 = threshold_image_2d(image1, i);
            std::vector<std::vector<bool>> thresh2 = threshold_image_2d(image2, j);

            std::vector< std::vector< bool > > binary_image_ij = and_two_binary_images(thresh1, thresh2);

            surface[i][j] = char_binary_image_2d(binary_image_ij);
        }
    }

    return surface;
}

//================================================

std::vector<std::vector<int>> surface_2d(const std::vector<std::vector<int>> &image1,
    const std::vector<std::vector<int>> &image2,
    int M,
    const std::vector<int> &vec_char)
{
    size_t num_rows = image1.size();
    size_t num_cols = image1[0].size();

    const std::vector<size_t> powers_of_two = {1, 2, 4, 8, 0,
                                               16, 32, 64, 128};
    // Check length of key
    if ( image1.size() != image2.size() || image1[0].size() != image2[0].size())
        throw std::runtime_error("Images dimensions must be equal.");

    // Output Euler char surface: numpy array
    std::vector<std::vector<int>> surface(M, std::vector<int>(M, 0));

    // Padded images
    std::vector<std::vector<int>> padded1 = pad(image1);
    std::vector<std::vector<int>> padded2 = pad(image2);

    // Loop over all pixel in original image
    for (size_t i = 1; i < num_rows+1; i++)
    {
        for (size_t j = 1; j < num_cols+1; j++)
        {
            int value2 = padded2[i][j];
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

            size_t num2 = num_after(padded2, indices);

            // Find elements in padded1[indices] greater than value of central_pixel
            int min_value = padded1[indices[4][0]][indices[4][1]];
            std::vector<int> euler_changes(M, 0);
            std::vector<int> sorted_values(10, M);
            for (int k = 0; k < 9; k++)
            {
                // Check value1 will produce some change
                // I am not checking if value1 was already iterated on
                int value1 = padded1[indices[k][0]][indices[k][1]];
                if (value1 < min_value)
                {
                    continue;
                }
                else
                {
                    sorted_values[k] = value1;
                }

                size_t num1 = 0;
                // Threshold padded1 on indices to obtain num1
                for (int s = 0; s < 9; s++)
                {
                    if ( padded1[indices[s][0]][indices[s][1]] <= value1 )
                        num1 |= powers_of_two[s];
                }

                size_t num = ( num1 & num2 );
                euler_changes[value1] = vec_char[num];
            }

            // Sort (and make unique) sorted_values and std::fill
            // euler_changes into surface
            std::sort(sorted_values.begin(), sorted_values.end());
            auto last = std::unique(sorted_values.begin(), sorted_values.end());
            sorted_values.erase(last, sorted_values.end());

            for (size_t t = 0; t < sorted_values.size()-1; t++)
            {
                size_t ind1 = static_cast<size_t>(sorted_values[t]);
                size_t ind2 = static_cast<size_t>(sorted_values[t+1]);
                for (size_t z = ind1; z < ind2; z++)
                    surface[value2][z] += euler_changes[ind1];
            }
        }
    }

    // Cumulative sum along columns
    std::vector<int> container = surface[0];
    for (size_t i = 1; i < M; i++)
    {
        for (size_t j = 0; j < M; j++)
        {
            surface[i][j] += container[j];
        }
        container = surface[i];
    }

    return surface;
}

