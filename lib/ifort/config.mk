# Voro++, a 3D cell-based Voronoi library
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 28th 2011

# This a common configuration file that includes definitions used by all
# the Makefiles.

# C++ compiler
CXX=icpc

ifeq ($(shell uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN),1)
     ARCH= -axSSSE3,SSE4.1,SSE4.2,AVX,CORE-AVX2 -mmacosx-version-min=10.12 -mdynamic-no-pic
else
     ARCH= -axSSE2,SSSE3,SSE4.1,SSE4.2,AVX
endif

# Flags for the C++ compiler
CFLAGS=-Wall -ansi -pedantic -O3 -no-prec-div -fp-model fast=2 -DVOROPP_VERBOSE=1 $(ARCH)

# Relative include and library paths for compilation of the examples
E_INC=-I../../src
E_LIB=-L../../src

# Installation directory
PREFIX=/usr/local

# Install command
INSTALL=install

# Flags for install command for executable
IFLAGS_EXEC=-m 0755

# Flags for install command for non-executable files
IFLAGS=-m 0644
