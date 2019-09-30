#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "libsurface.hpp"

namespace py = pybind11;

PYBIND11_MODULE(surface, m) {
    m.doc() = "Euler characteristic surfaces cpp bindings.";

    m.def("naive_surface_2d", &naive_surface_2d,
          "Euler char surface of 2d images with naive algorithm.",
          py::arg("image1"),
          py::arg("image2"),
          py::arg("max_intensity"));
    m.def("surface_2d", &surface_2d,
          "Euler char surface of 2d images.",
          py::arg("image1"),
          py::arg("image2"),
          py::arg("max_intensity"),
          py::arg("vector_euler_changes"));
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
