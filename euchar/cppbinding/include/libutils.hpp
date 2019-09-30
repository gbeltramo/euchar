#ifndef EUCHAR_LIBUTILS_H
#define EUCHAR_LIBUTILS_H

int sum_bool_matrix(std::vector<std::vector<bool>> matrix);

std::vector<std::vector<bool>> threshold_image_2d(const std::vector<std::vector<int>> &image, int value);

std::vector<std::vector<bool>> and_two_binary_images(const std::vector<std::vector<bool>> &binary1, const std::vector<std::vector<bool>> &binary2);

int char_binary_image_2d(std::vector<std::vector<bool>> input);

std::vector<std::vector<bool>> neigh_2d_from_number(uint8_t num);

uint8_t number_from_neigh_2d(std::vector<std::vector<bool>> neigh);

std::vector<short> vector_of_euler_changes_2d();

int sum_bool_matrix_3d(std::vector<std::vector<std::vector<bool>>> matrix);

int char_binary_image_3d(std::vector<std::vector<std::vector<bool>>> input);

std::vector<std::vector<std::vector<bool>>> neigh_3d_from_number(size_t num);

size_t number_from_neigh_3d(std::vector<std::vector<std::vector<bool>>> neigh);

std::vector<short> vector_of_euler_changes_3d(size_t offset=0, size_t max_value=6710886);

#endif