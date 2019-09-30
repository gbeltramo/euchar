#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "libcurve.hpp"

namespace py = pybind11;

PYBIND11_MODULE(curve, m) {
    m.doc() = "Euler characteristic curves cpp bindings.";

    m.def("pad", &pad,
	  "Padded image.");
    m.def("num_after", &num_after,
	  "Number of euler change after adding central pixel.");
    m.def("naive_image_2d", &naive_image_2d,
	  "Euler char curve with nested loop.");
    m.def("image_2d", &image_2d,
	  "Euler char curve of 2d image.");
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
