cmake_minimum_required(VERSION 3.5)
project(module_cpp_binding VERSION 0.1
                           LANGUAGES CXX)

find_package(pybind11)

# UTILS

pybind11_add_module(libutils STATIC src/libutils.cpp)
target_compile_features(libutils PUBLIC cxx_std_11)

pybind11_add_module(utils utils.cpp include/libutils.hpp)
target_include_directories(utils PUBLIC include/)
target_link_libraries(utils PUBLIC libutils)

# CURVE
pybind11_add_module(libcurve STATIC src/libcurve.cpp include/libutils.hpp)
target_compile_features(libcurve PUBLIC cxx_std_11)
target_include_directories(libcurve PUBLIC include/)
target_link_libraries(libcurve PUBLIC libutils)

pybind11_add_module(curve curve.cpp include/libcurve.hpp)
target_include_directories(curve PUBLIC include/)
target_link_libraries(curve PRIVATE libcurve)

# SURFACE
pybind11_add_module(libsurface STATIC src/libsurface.cpp
include/libutils.hpp)
target_compile_features(libsurface PUBLIC cxx_std_11)
target_include_directories(libsurface PUBLIC include/)
target_link_libraries(libsurface PUBLIC libutils)

pybind11_add_module(surface surface.cpp include/libsurface.hpp)
target_include_directories(surface PUBLIC include/)
target_link_libraries(surface PRIVATE libsurface)
