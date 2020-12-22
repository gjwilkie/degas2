# Notes for converting degas2 build system:
# - Suffix fortran90 files with .f90 

# Must init the autoconf setup
# The first parameter is project name
# second is version number
# third is bug report address

AC_INIT([degas2], [1.0])

# Safety checks in case user overwritten --srcdir
AC_CONFIG_SRCDIR([randomtest.web])
AC_CONFIG_MACRO_DIR([m4])

# Need a c compiler for some parts
AC_PROG_CC

AX_PROG_F77_MPI([],[],[AC_MSG_WARN([Could not find MPI, disabling.]);AC_PROG_F77])
#AX_PROG_F77_MPI()

AC_FC_FIXEDFORM([],[AC_MSG_ERROR([Fortran compiler does not support fixed form, which is required by fweb output.])])

# Argument mismatch is used a lot in degas2. Newer versions of gfortran interpret this as an error
# unless the following flag is defined. Check if it's needed.
CHECK_COMPILE_FLAG(Fortran,[-fallow-argument-mismatch],FFLAGS,mismatch_flag)

#FFLAGS+= -fPIC -I../silo-4.10.2/include triangle.o 

# Include std math library
AC_CHECK_LIB([m],[cos])

#AC_SEARCH_LIBS([silo],[dbcreate],[],[AC_MSG_ERROR([Could not find SILO library, required by several tools including definegeometry2d.])

# Include netcdf library
AC_CHECK_LIB(netcdff,ncopn_,[],[AC_MSG_ERROR([unable to find netcdff library; required])])
PKG_CHECK_MODULES([LIBNETCDF],[netcdf-fortran])

# Other packages:
# SILO: needed for plotgeom, geomtesta, and ucd_plot. Distributed.
# Triangle: as git submodule
# HDF5 is a possibility; influence which SILO option to use. Checking for zlib, etc. is complicated (see configure.ac in silo package). Leaving this out for now and just setting GRAPH_FILE=SILO

# FWEB flags to define:
# MPI=yes/no
# COMPILER
# SYSTEM
# GRAPH_FILE
# DEBUG
# MASTER_SAMPLE


#AX_MPI([mpi_found=true],[mpi_found=false])

# Init automake, and specify this program use relaxed structures.
# i.e. this program doesn't follow the gnu coding standards, and doesn't have
# ChangeLog, COPYING, AUTHORS, INSTALL, README etc. files.
#AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_INIT_AUTOMAKE([-Wall foreign])


# Tells automake to create a Makefile
# See https://www.gnu.org/software/automake/manual/html_node/Requirements.html
AC_CONFIG_FILES([Makefile])

# Generate the output
AC_OUTPUT
