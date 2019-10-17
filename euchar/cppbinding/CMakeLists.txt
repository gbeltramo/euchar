cmake_minimum_required(VERSION 3.8)
project(module_cpp_binding VERSION 0.1
                           LANGUAGES CXX)

find_package(pybind11)

pybind11_add_module(utils MODULE utils.cpp)
target_compile_features(utils PUBLIC cxx_std_11)

pybind11_add_module(curve MODULE curve.cpp)
target_compile_features(curve PUBLIC cxx_std_11)

pybind11_add_module(surface MODULE surface.cpp)
target_compile_features(surface PUBLIC cxx_std_11)
