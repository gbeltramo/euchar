#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include "libsurface.hpp"

namespace py = pybind11;


PYBIND11_MODULE(surface, m) {
    m.doc() = "Euler characteristic surfaces cpp bindings.";

    m.def("naive_images_2d", &naive_images_2d,
          "Euler char surface of 2d images with naive algorithm.");
    m.def("images_2d", &images_2d,
          "Euler char surface of 2d images");
    m.def("naive_images_3d", &naive_images_3d);
    m.def("images_3d", &images_3d);
    m.def("bifiltration", &bifiltration,
          "Euler char surface of 2 filtrations",
          py::arg("simplices"),
          py::arg("parametrization1"), py::arg("parametrization2"),
          py::arg("bins1"),py::arg("bins2"));
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
