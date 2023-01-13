# libpebble

Calculate higher order equilibrium classes of 3D scanned pebbles.

## Directory structure

The source code is organized using the following directory structure:
- `build`: binary output directory, see the Compilation section below for more information on how to build the binaries
  - `conanenv`: contains the files required to set up the build environment
- `example`: contains examples of using the library
  - `pebble-cli.cpp`: this example takes the 3D scan and the primary class of a pebble and calculates its secondary and tertiary classes along with numerous other shape descriptors. Run the example with `pebble-cli -h` for the complete list of input and output values.
  - `graph2tikz.cpp`: this applications takes the alphanumerical encoding of a Reeb, Morse-Smale or master graph output by `pebble-cli` and draws its planar embedding in LaTeX format using the TikZ library
- `include/pebble`: the source files of the header-only `libpebble` library
- `test`: unit tests for some of the functions in the library

## Compilation

You need an up-to-date C++ compiler and Python 3 installation on your system in order to build the project. Every other build tool and library dependency will be installed using [Conan](https://docs.conan.io). The code snippet below goes through the following steps:
1. Create and activate a Python virtual environment.
2. Install Conan in the virtual environment.
3. Create and activate a Conan virtual environment. This makes sure that CMake, Ninja and Doxygen are available to build the library.
4. Install all the required libraries (boost, CGAL, bliss, GoogleTest) using Conan.
5. Configure the project using CMake then build it with Ninja.

```
cd build
python3 -m venv pyenv
source pyenv/bin/activate
pip install conan
cd conanenv
conan install .
source activate.sh
cd ..
conan install ..
cmake .. -GNinja
ninja
```

Among the generated directories are the following:
- `example`: directory containing the `pebble-cli` and `graph2tikz` executables
- `html`: the documentation of the library functions
- `test`: test cases for the library
