#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "libutils.hpp"

namespace py = pybind11;
  
//================================================

PYBIND11_MODULE(utils, m) {
    m.doc() = "";
    m.def("sum_bool_matrix", &sum_bool_matrix,
	  "Number of True values in bool_matrix",
	  py::arg("bool_matrix")
	  );
    m.def("char_binary_image_2d", &char_binary_image_2d);
    m.def("vector_of_euler_changes_2d", &vector_of_euler_changes_2d,
	  "Vector of all possible Euler characteristic changes for 2D images");
    m.def("sum_bool_matrix_3d", &sum_bool_matrix_3d,
          "docs sum bool matrix 3d",
          py::arg("bool_matrix_3d"));
    m.def("char_binary_image_3d", &char_binary_image_3d);
    m.def("neigh_3d_from_number", &neigh_3d_from_number);
    m.def("number_from_neigh_3d", &number_from_neigh_3d);
    m.def("vector_of_euler_changes_3d", &vector_of_euler_changes_3d,
	  "Vector of all possible Euler characteristic changes for 3D images");
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
