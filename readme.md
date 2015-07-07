# mpm-2d

An implementation of the material point method (MPM) in 2D, with particular focus on granular materials. It is not feature complete for general use, as its primary purpose is research code to test material models in situations which are difficult for current modeling techniques.

`mpm-2d` consists of several parts:

 * MPM library
 * Material models
 * Driver program for a reference MPM implementation
 * Visualization program for the output of the reference implementation
 * Examples for the reference program as well as supporting scripts


## Building

`mpm-2d` uses [CMake](http://www.cmake.org/) from Kitware as its build tool.
You will need at least CMake 2.8.12. You can build from source by running the following commands.

```
mkdir build
cd build
cmake ..
make
```

After these steps, you can do a `make install`. This will move binaries and examples into the appropriate place. The default is under the build directory in the `installed` dir. You can set the build dir by replacing the `cmake ..` command with `cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..` but leave the trailing space off. The binaries will be located under the `bin` directory, while data files such as the examples and material models will be under the `mpm` directory. The driver program is `mpm_2d` and the visualization program is `mpm_viz`.

## Running

By default, the driver looks for a configuration file `simulation.cfg` in the current working directory. This file specifies initial conditions, geometry, output parameters, and the material model (see the example configurations for more details). you can override certain options on the command line, which are shown when calling `mpm_2d -h`. The most important of these is the specification of the configuration file given by `mpm_2d -c CFGFILE`, however the command line arguments (if specified) will override the file paramters.

## Dependencies

Several libraries are used in creation of the driver program, as well as some of the material models. CMake will complain if you do not have the required libraries on your system, but lack of the optional libraries will mean that the visualization program will not be built.

 * [CXSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) - Required, or equivalent such as suitesparse (same page).
 * [libconfuse](http://www.nongnu.org/confuse/) - Required.
 * [FTGL](http://sourceforge.net/projects/ftgl/) - Optional.
 * [FreeType](http://www.freetype.org/index.html) - Optional.
 * [OpenGL](https://www.opengl.org/) - Optional.
 * [png](http://www.libpng.org/pub/png/) - Optional.
 * [SDL](https://www.libsdl.org/) - Optional.
