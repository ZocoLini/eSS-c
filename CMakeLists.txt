cmake_minimum_required(VERSION 3.10)

project(ess-c LANGUAGES C Fortran)

set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort)

set(CMAKE_C_FLAGS "-O3 -ipo -xHost -DINTEL -diag-disable=10441")
set(CMAKE_Fortran_FLAGS "-O3 -ipo -xHost -fpp -DEXPORT -DINTEL")

# set(CLIBS -cxxlib -lrt -lhdf5 -lxml2 -limf -lifcore -lm -MMD -lz)
# set(FLIBS -cxxlib)
# set(LIBS -L${PROJECT_SOURCE_DIR}/lib/BLAS -lblas)
# set(INC -I${PROJECT_SOURCE_DIR}/lib/misqp/gnu)

add_subdirectory(sacess)
