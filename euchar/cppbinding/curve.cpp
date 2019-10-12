#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "libcurve.hpp"

namespace py = pybind11;

PYBIND11_MODULE(curve, m) {
    m.doc() = "Euler characteristic curves cpp bindings.";

    // images functions
    m.def("naive_image_2d", &naive_image_2d,
	  "Euler char curve of 2d image with nested loops on pixels.");
    m.def("image_2d", &image_2d,
	  "Euler char curve of 2d image.");
    m.def("naive_image_3d", &naive_image_3d,
          "Euler char curve of 3d image with nested loops on pixels.");
    m.def("image_3d", &image_3d,
          "Euler char curve of 3d image.");
    
    // filtrations functions
    m.def("filtration", &filtration);
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
