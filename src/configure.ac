# Notes for converting degas2 build system:
# - Suffix fortran90 files with .f90 

# Must init the autoconf setup
# The first parameter is project name
# second is version number
# third is bug report address
AC_INIT([degas2], [1.0])

# Safety checks in case user overwritten --srcdir
AC_CONFIG_SRCDIR([randomtest.web])

# Store the auxiliary build tools (e.g., install-sh, config.sub, config.guess)
# in this dir (build-aux)
AC_CONFIG_AUX_DIR([build-aux])

AC_PROG_F77

AC_FC_FIXEDFORM([],[AC_MSG_ERROR([Fortran compiler does not support fixed form, which is required by fweb output.])])
AC_CHECK_LIB(netcdff,ncopn_,[],[AC_MSG_ERROR([unable to find netcdff library])])
PKG_CHECK_MODULES([LIBNETCDF],[netcdf-fortran])

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
