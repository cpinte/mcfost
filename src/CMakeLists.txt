cmake_minimum_required(VERSION 3.12)
project(MCFOST_CPP_Interface LANGUAGES CXX Fortran)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)

# Find OpenMP for parallel execution
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(mcfost_interface PUBLIC OpenMP::OpenMP_CXX)
endif()
if(OpenMP_Fortran_FOUND)
    target_link_libraries(mcfost_interface PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Find MPI if available (MCFOST may use MPI)
find_package(MPI)
if(MPI_FOUND)
    add_definitions(-DUSE_MPI)
endif()

# Compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -cpp -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=bounds -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -fpp -extend-source")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -check bounds -check all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -unroll")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -pedantic")
    set(CMAKE
