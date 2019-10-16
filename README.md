# euchar

Compute Euler Characteristic curves and surfaces of image and point data.

## Installation

For the package to install correctly the prerequisites are

- CMake version >= 3.5
- pybind11 version >= 2.2

This module makes use of features provided by the following Python modules

- numpy
- matplotlib, used in the display module of the package
- scipy, used to obtain Delaunay triangulations of 2D and 3d finite point sets
- scikit-learn, for nearest neighbours algorithm to obtain an estimate of the inverse of the local density in a finite set of points


## Examples

See the `/notebooks/` directory of the github page for usage examples.


## Documentation

The documentation for this project is available at https://gbeltramo.github.io/euchar/