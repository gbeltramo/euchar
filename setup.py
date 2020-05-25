import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = sourcedir


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.sourcedir)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', os.path.abspath(ext.sourcedir)] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

#===========================================================

with open("README.md", "r") as f:
    long_description = f.read()

utils = CMakeExtension(name="utils",
                       sourcedir="euchar/cppbinding/")

curve = CMakeExtension(name="curve",
                       sourcedir="euchar/cppbinding/")

surface = CMakeExtension(name="surface",
                         sourcedir="euchar/cppbinding/")
setup(
    name='euchar',
    version='0.1.6',
    packages=find_packages("."),
    author="Gabriele Beltramo",
    author_email="gabri.beltramo@gmail.com",
    description="Euler characteristic of images and finite point sets.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gbeltramo/euchar",
    license='MIT',
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'scikit-learn',
        'pybind11'
    ],
    ext_modules=[utils, curve, surface],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    python_requires='>=2.7,!=3.1,!=3.2,!=3.3,!=3.4',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='topology data analysis, euler characteristic, euler characteristic curve, euler characteristic surface'
)
