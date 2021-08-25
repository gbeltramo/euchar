# euchar

Compute Euler Characteristic curves and surfaces of image and point data.

## Installation

The `setuptools`, `numpy`, and `scipy` Python packages are
prerequisites for using this package.

You'll also need `pybind11`, `cmake` and a `C++` compiler in order to install this
package. It is recommended to install both `pybind11` and `cmake` via `conda`.

```
>>> conda install -c anaconda cmake
>>> conda install -c conda-forge pybind11
```

Finally run the following command to build and install this Python package

```
>>> pip install euchar
```

If the command above fails to build this package, then you might want to clone
this repository to a directory `euchar/` on your computer and run

```
>>> cd /<path>/<to>/<cloned>/<repo>/euchar/
>>> python setup.py install
```

which outputs more information on the errors causing the build to fail.

**Windows.** After installing `conda`, run the above commands within an
`Anaconda prompt`. For the C++ compiler install
<a href="https://visualstudio.microsoft.com/vs/">Visual Studio community</a>.

## Documentation

The documentation for this project is available at https://gbeltramo.github.io/euchar/,
where usage examples are provided for computing both Euler characteristic
curves and Euler characteristic surfaces.
