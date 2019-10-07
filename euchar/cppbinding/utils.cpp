#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "libutils.hpp"

namespace py = pybind11;
  
//================================================

PYBIND11_MODULE(utils, m) {
    m.doc() = "";

    // 2D functions
    m.def("sum_bool_2d", &sum_bool_2d,
	  "Number of True values in bool_matrix",
	  py::arg("bool_matrix")
	  );
    m.def("pad_2d", &pad_2d);
    m.def("threshold_image_2d", &threshold_image_2d);
    m.def("elementwise_AND_2d", &elementwise_AND_2d);
    m.def("neigh_pixel_2d", &neigh_pixel_2d);
    m.def("binary_neigh_pixel_2d", &binary_neigh_pixel_2d);
    m.def("char_binary_image_2d", &char_binary_image_2d);
    m.def("neigh_2d_from_number", &neigh_2d_from_number);
    m.def("number_from_neigh_2d", &number_from_neigh_2d);
    m.def("vector_of_euler_changes_2d", &vector_of_euler_changes_2d,
	  "Vector of all possible Euler characteristic changes for 2D images");

    // 3D functions
    m.def("sum_bool_3d", &sum_bool_3d,
          "docs sum bool matrix 3d",
          py::arg("bool_matrix_3d"));
    m.def("pad_3d", &pad_3d,
          "Returns padded 3d image",
          py::arg("image_3d"), py::arg("max_intensity"));
    m.def("threshold_image_3d", &threshold_image_3d,
          "Threshold of 3d image",
          py::arg("image_3d"), py::arg("threshold"));
    m.def("elementwise_AND_3d", &elementwise_AND_3d);
    m.def("neigh_voxel_3d", &neigh_voxel_3d);
    m.def("binary_neigh_voxel_3d", &binary_neigh_voxel_3d);
    m.def("char_binary_image_3d", &char_binary_image_3d);
    m.def("neigh_3d_from_number", &neigh_3d_from_number);
    m.def("number_from_neigh_3d", &number_from_neigh_3d);
    m.def("vector_of_euler_changes_3d", &vector_of_euler_changes_3d,
	  "Vector of all possible Euler characteristic changes for 3D images");

    // points utils
    m.def("dim_simplex", &dim_simplex);
    m.def("filter_parametrization_on", &filter_parametrization_on);
    
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
