## DEGAS2 fork
This is a fork of DEGAS 2 with some installation and running scripts (which as of now are incomplete). 

### Installation/build
1. make a build directory (e.g. inside of the degas2 repo.),
    ```
    mkdir build
    cd build
    ```
2. Run CMake,
    ```
    cmake ..
    ```
    *Optional*: If you wish to enable MPI (which currently requires flighttest to be run with mpirun), use the following CMake option:
    ```
    cmake -DUSE_MPI=ON ..
    ```
    *Optional*: If you wish to enable debug mode use: `cmake -DCMAKE_BUILD_TYPE=Debug ..`

3. Make and test `datasetup` program, 
    From the build directory,
    ```
    make datasetup
    cp ../data/degas2.in .
    bin/datasetup
    ```
4. Make the remaining programs, 
    From the build directory,
    ```
    make boxgen problemsetup definegeometry2d defineback tallysetup flighttest ucd_plot outputbrowser
    ```

### Running simulations

The degas2.in and tally.in files in degas2/data/ should be copied to each of your working directories. 

Executables are created in the build/bin and can be run from an arbitrary working directory if this directory is in your PATH environment variable. You may also find it handy to include degas2/scripts in your PYTHONPATH environment variable.

### Notes for different platforms

If building on MacOS, FWEB won't compile with the default "gcc" (clang wrapper). Install the GNU version of gcc from macports or homebrew and replace the above cmake command with:
   cmake -DCMAKE_C_COMPILER=gcc-11 ..
(or whichever executable and its path is GNU-gcc). This additional option only needs to be run the first time you configure on MacOS.

*Update Fall 2026*: With newer versions of GNU-gcc one might need to specify other arguments for cmake. 
I was able to get a succesful build on macOS 26.7 from scratch using the following brew-install commands and cmake options,
```sh
conda deactivate # in case a conda env is active.
brew install cmake
brew install gcc
brew install open-mpi
cmake -DCMAKE_C_COMPILER=gcc-15 \
      -DCMAKE_CXX_COMPILER=g++-15 \
      -DUSE_MPI=ON \
      -DCMAKE_C_FLAGS="-Dstat64=stat -std=gnu17" \
      -DMPI_C_COMPILER=/opt/homebrew/bin/mpicc \
      -DMPI_CXX_COMPILER=/opt/homebrew/bin/mpicxx \
      -DMPI_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
      ..
```

On the PPPL cluster, ensure the following modules are loaded: gcc, openmpi, szip, hdf5-parallel, netcdf-c, netcdf-fortran, silo, cmake. And define the environment variable "export NetCDF_ROOT=$NETCDF_C_HOME:$NETCDF_FORTRAN_HOME".

On the `omega` GA cluster, make sure `degas2_build` module is loaded. 

### Other

All relevant documentation on compiling and running DEGAS 2 is 
contained in the manual (Doc/degas2_all.pdf). See the section there entitled "Compiling DEGAS 2".

DEGAS 2  usage is currently governed by a licensing procedure.  Interested
users should complete and submit the form:

https://docs.google.com/forms/d/e/1FAIpQLSdiPXgLqcmZfiX5sEqf6TOVHd8igK-0v-gBHMFa-mBsPn6QIQ/viewform


