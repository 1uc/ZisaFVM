# ZisaFVM
[![Build Status](https://github.com/1uc/ZisaFVM/actions/workflows/basic_integrity_checks.yml/badge.svg?branch=main)](https://github.com/1uc/ZisaFVM/actions)
[![Docs Status](https://github.com/1uc/ZisaFVM/actions/workflows/publish_docs.yml/badge.svg?branch=main)](https://1uc.github.io/ZisaFVM)

ZisaFVM is an unstructured, high-order finite volume solver used to research
a well-balanced modification of the base scheme.

## Quickstart
Start by cloning the repository

    $ git clone https://github.com/1uc/ZisaFVM.git

and change into the newly created directory. Then proceed to install the
dependencies:

    $ bin/install_dir.sh COMPILER DIRECTORY DEPENDENCY_FLAGS

they will be placed into a subdirectory of `DIRECTORY` and print
part of the CMake command needed to include the dependencies. `COMPILER` must
be replaced with the compiler you want to use. The available `DEPENDENCY_FLAGS`
are

  * `--zisa_has_mpi={0,1}` to request MPI.
  * `--zisa_has_netcdf={0,1}` ZisaFVM uses HDF5 for I/O. However, lower level
  libraries meanwhile support NetCDF. You can build the dependencies with NetCDF
  if you need it yourself. Otherwise, there is no need to request NetCDF.

If this worked continue by running the `cmake` command and compiling the
library. Take a look at the [project specific flags] for CMake if you want to
modify something. If this didn't work, it's not going to be a quick start.
Please read [Dependencies] and then [Building].

[project specific flags]: https://1uc.github.io/ZisaFVM/md_building.html#cmake_flags
[Dependencies]: https://1uc.github.io/ZisaFVM/md_dependencies.html
[Building]: https://1uc.github.io/ZisaFVM/md_building.html
