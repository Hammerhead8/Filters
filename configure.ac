#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Filters], [1.0], [josh.cates@protonmail.com])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([filters.cpp])
AC_CONFIG_HEADERS([config.h])

# Set the language to C++
AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Initialize libtool
LT_INIT([shared])

# Check for the math library
AC_CHECK_LIB(m, cos)

# Check for the math and vector header files
AC_CHECK_HEADERS([cmath vector])

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt cos asinh sinh cosh])

# Create the makefile
AC_CONFIG_FILES([Makefile])

AC_OUTPUT
