# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

project(ess-c LANGUAGES C Fortran)

set(BUILDDIR .)
set(LIBRARY ${BUILDDIR}/lib)
set(AMIGO_PATH ${LIBRARY}/libAMIGO)

set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort)

set(CMAKE_C_FLAGS "-O3 -ipo -xHost -DINTEL")
set(CMAKE_Fortran_FLAGS "-O3 -ipo -xHost -fpp -DEXPORT -DINTEL")

set(CLIBS -cxxlib -lrt -lhdf5  -lxml2 -limf -lifcore  -lm -MMD -lz)
set(FLIBS -cxxlib)
set(LIBS -L${LIBRARY}/BLAS -lblas)
set(INC -I${LIBRARY}/misqp/gnu)

file(GLOB SRC_C
        ${BUILDDIR}/benchmarks/bbob/*.c
        ${BUILDDIR}/benchmarks/customized/*.c
        ${BUILDDIR}/benchmarks/systemsBiology/*.c
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B1/*.c
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B2/*.c
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B3/*.c
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B4/*.c
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B5/*.c
        ${BUILDDIR}/src/input_module/*.c
        ${BUILDDIR}/src/output/*.c
        ${BUILDDIR}/src/structure_data/*.c
        ${BUILDDIR}/src/method_module/*.c
        ${BUILDDIR}/src/method_module/DE/*.c
        ${BUILDDIR}/src/method_module/eSS/*.c
        ${BUILDDIR}/src/*.c
)

file(GLOB SRC_FORTRAN
        ${BUILDDIR}/src/method_module_fortran/eSS/*.f90
        ${BUILDDIR}/src/method_module_fortran/eSS/localsolvers/*.f90
)

file(GLOB SRC_AMIGO
        ${AMIGO_PATH}/src/src_cvodes/*.c
        ${AMIGO_PATH}/src/src_amigo/*.c
        ${AMIGO_PATH}/src/src_de/*.c
        ${AMIGO_PATH}/src/src_fullC/*.c
        ${AMIGO_PATH}/src/src_julia/*.c
        ${AMIGO_PATH}/src/src_SRES/*.c
)

include_directories(
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B1/
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B6/fly/
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B6/util/
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B6/
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B5/
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B4/
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B3/
        ${BUILDDIR}/benchmarks/systemsBiology/BioPredyn/B2/
        ${BUILDDIR}/benchmarks/systemsBiology/others/circadian/
        ${BUILDDIR}/benchmarks/systemsBiology/others/Nfkb/
        ${BUILDDIR}/benchmarks/systemsBiology/others/3-step_pathway/
        ${AMIGO_PATH}/include/include_cvodes/
        ${AMIGO_PATH}/include/include_amigo/
        ${AMIGO_PATH}/include/include_de/
        ${AMIGO_PATH}/include/include_SRES/
        ${BUILDDIR}/include/
        ${BUILDDIR}/include/input/
        ${BUILDDIR}/include/method_module/
        ${BUILDDIR}/include/method_module/eSS/
        ${BUILDDIR}/include/method_module/DE/
        ${BUILDDIR}/include/bbob/
        ${BUILDDIR}/include/structure_data/
        ${BUILDDIR}/include/output/
        ${BUILDDIR}/include/error/
        ${BUILDDIR}/include/method_module_fortran/
        ${BUILDDIR}/src/method_module_fortran/eSS/
        ${BUILDDIR}/src/method_module_fortran/eSS/localsolvers
        ${BUILDDIR}/src/method_module_fortran/eSS/localsolvers/dhc
        ${BUILDDIR}/src/method_module_fortran/eSS/localsolvers/misqp
        ${BUILDDIR}/src/method_module_fortran/eSS/localsolvers/nl2sol
        ${LIBRARY}/gsl-1.14/include
)

link_directories(${LIBRARY}/gsl-1.14/lib)
link_directories(${LIBRARY}/libxml2/lib)
link_directories(${LIBRARY}/hdf5-1.8.12/lib)
link_directories(${AMIGO_PATH}/lib)

add_executable(ess-c
        ${SRC_C}
        ${SRC_FORTRAN}
        ${SRC_AMIGO}
)

add_library(libAMIGO STATIC IMPORTED)
set_target_properties(libAMIGO PROPERTIES IMPORTED_LOCATION ${AMIGO_PATH}/lib/libAMIGO.a)

target_link_libraries(ess-c
        ${LIBS}
        ${CLIBS}
        ${FLIBS}
        libAMIGO
)
