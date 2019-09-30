#ifndef EUCHAR_LIBSURFACE_H
#define EUCHAR_LIBSURFACE_H

std::vector<std::vector<int>> naive_surface_2d(const std::vector<std::vector<int>> &image1, const std::vector<std::vector<int>> &image2, int M);

std::vector<std::vector<int>> surface_2d(const std::vector<std::vector<int>> &image1, const std::vector<std::vector<int>> &image2, int M, const std::vector<int> &vec_char);

#endif