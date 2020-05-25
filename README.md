# euchar

Compute Euler Characteristic curves and surfaces of image and point data.

## Installation

For the package to install correctly the prerequisites are

- a `C++` compiler. Windows users need to install a recent version of
[Visual Studio](https://visualstudio.microsoft.com/vs/)
- `cmake` version >= 3.8
- `pybind11` version >= 2.2

It is recommended to install both `cmake` and `pybind11` via the `conda` command,
which can be obtained by installing either [Anaconda](https://www.anaconda.com/)
or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

With `conda` installed run the following commands in a terminal window.

```
$ conda install -c anaconda cmake
$ conda install -c conda-forge pybind11
```

This module makes use of features provided by the following Python modules

- numpy
- matplotlib, used in the display module of the package
- scipy, used to obtain Delaunay triangulations of 2D and 3d finite point sets
- scikit-learn, for nearest neighbours algorithm to obtain an estimate of the inverse of the local density in a finite set of points

## Documentation

The documentation for this project is available at https://gbeltramo.github.io/euchar/,
where usage examples are provided for computing both Euler characteristic
curves and Euler characteristic surfaces.
