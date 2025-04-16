# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

project(ess-c LANGUAGES C Fortran)

set(BUILDDIR .)
set(LIBRARY ${PROJECT_SOURCE_DIR}/lib)
set(AMIGO_PATH ${PROJECT_SOURCE_DIR}/lib/libAMIGO)

set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort)

set(CMAKE_C_FLAGS "-O3 -ipo -xHost -DINTEL -diag-disable=10441")
set(CMAKE_Fortran_FLAGS "-O3 -ipo -xHost -fpp -DEXPORT -DINTEL")

set(CLIBS -cxxlib -lrt -lhdf5 -lxml2 -limf -lifcore -lm -MMD -lz)
set(FLIBS -cxxlib)
set(LIBS -L${LIBRARY}/BLAS -lblas)
set(INC -I${LIBRARY}/misqp/gnu)

file(GLOB SRC_C
        ${PROJECT_SOURCE_DIR}/benchmarks/bbob/*.c
        ${PROJECT_SOURCE_DIR}/benchmarks/customized/*.c
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/*.c
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B1/*.c
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B2/*.c
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B3/*.c
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B4/*.c
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B5/*.c
        ${PROJECT_SOURCE_DIR}/src/input_module/*.c
        ${PROJECT_SOURCE_DIR}/src/output/*.c
        ${PROJECT_SOURCE_DIR}/src/structure_data/*.c
        ${PROJECT_SOURCE_DIR}/src/method_module/*.c
        ${PROJECT_SOURCE_DIR}/src/method_module/DE/*.c
        ${PROJECT_SOURCE_DIR}/src/method_module/eSS/*.c
        ${PROJECT_SOURCE_DIR}/src/*.c
)

file(GLOB SRC_FORTRAN
        ${PROJECT_SOURCE_DIR}/src/method_module_fortran/eSS/*.f90
        ${PROJECT_SOURCE_DIR}/src/method_module_fortran/eSS/localsolvers/*.f90
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
        ${PROJECT_SOURCE_DIR}/benchmarks/
        ${PROJECT_SOURCE_DIR}/benchmarks/bbob/
        ${PROJECT_SOURCE_DIR}/benchmarks/customized/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B1/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B6/fly/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B6/util/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B6/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B5/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B4/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B3/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/BioPredyn/B2/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/others/circadian/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/others/Nfkb/
        ${PROJECT_SOURCE_DIR}/benchmarks/systemsBiology/others/3-step_pathway/

        ${PROJECT_SOURCE_DIR}/include/
        ${PROJECT_SOURCE_DIR}/include/input/
        ${PROJECT_SOURCE_DIR}/include/method_module/
        ${PROJECT_SOURCE_DIR}/include/method_module/eSS/
        ${PROJECT_SOURCE_DIR}/include/method_module/DE/
        ${PROJECT_SOURCE_DIR}/include/bbob/
        ${PROJECT_SOURCE_DIR}/include/structure_data/
        ${PROJECT_SOURCE_DIR}/include/output/
        ${PROJECT_SOURCE_DIR}/include/error/
        ${PROJECT_SOURCE_DIR}/include/method_module_fortran/
        ${PROJECT_SOURCE_DIR}/src/method_module_fortran/eSS/
        ${PROJECT_SOURCE_DIR}/src/method_module_fortran/eSS/localsolvers
        ${PROJECT_SOURCE_DIR}/src/method_module_fortran/eSS/localsolvers/dhc
        ${PROJECT_SOURCE_DIR}/src/method_module_fortran/eSS/localsolvers/misqp
        ${PROJECT_SOURCE_DIR}/src/method_module_fortran/eSS/localsolvers/nl2sol

        ${AMIGO_PATH}/include/
        ${AMIGO_PATH}/include/include_cvodes/
        ${AMIGO_PATH}/include/include_amigo/
        ${AMIGO_PATH}/include/include_de/
        ${AMIGO_PATH}/include/include_SRES/
        ${PROJECT_SOURCE_DIR}/lib/gsl-1.14/include
        ${PROJECT_SOURCE_DIR}/lib/libxml2/include
        ${PROJECT_SOURCE_DIR}/lib/hdf5-1.8.12/include
)

link_directories(${PROJECT_SOURCE_DIR}/lib/gsl-1.14/lib)
link_directories(${PROJECT_SOURCE_DIR}/lib/libxml2/lib)
link_directories(${PROJECT_SOURCE_DIR}/lib/hdf5-1.8.12/lib)
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
