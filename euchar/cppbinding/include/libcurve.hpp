#ifndef EUCHAR_LIBCURVE_H
#define EUCHAR_LIBCURVE_H

std::vector<std::vector<int>> pad(const std::vector<std::vector<int>> &image);

size_t num_after(const std::vector<std::vector<int>> &padded, const std::vector<std::vector<size_t>> &indices);

std::vector<int> naive_image_2d(std::vector<std::vector<int>> image, int M);

std::vector<int> image_2d(const std::vector<std::vector<int>> & image, int M, const std::vector< int > &vec_char);

#endif