file(GLOB SRC_C
        src/input_module/*.c
        src/output/*.c
        src/structure_data/*.c
        src/method_module/*.c
        src/method_module/DE/*.c
        src/method_module/eSS/*.c
        src/*.c

        benchmarks/bbob/*.c
        benchmarks/customized/*.c
        # benchmarks/systemsBiology/*.c
        # benchmarks/systemsBiology/BioPredyn/B1/*.c
        # benchmarks/systemsBiology/BioPredyn/B2/*.c
        # benchmarks/systemsBiology/BioPredyn/B3/*.c
        # benchmarks/systemsBiology/BioPredyn/B4/*.c
        # benchmarks/systemsBiology/BioPredyn/B5/*.c

        libAMIGO/src/src_cvodes/*.c
        libAMIGO/src/src_amigo/*.c
        libAMIGO/src/src_de/*.c
        libAMIGO/src/src_fullC/*.c
        libAMIGO/src/src_julia/*.c
        libAMIGO/src/src_SRES/*.c
)

file(GLOB SRC_FORTRAN
        src/method_module_fortran/eSS/*.f90
        src/method_module_fortran/eSS/localsolvers/*.f90
        src/method_module_fortran/eSS/localsolvers/dhc/*.f90
        src/method_module_fortran/eSS/localsolvers/misqp/*.f90
        src/method_module_fortran/eSS/localsolvers/nl2sol/*.f90
)

add_library(sacess
        ${SRC_C}
        ${SRC_FORTRAN}
)

include_directories(sacess
        include
        benchmarks
        libAMIGO/include
        libAMIGO/include/include_cvodes

        include/method_module_fortran

)

include_directories(
        ${PROJECT_SOURCE_DIR}/lib/gsl-1.14/include
        ${PROJECT_SOURCE_DIR}/lib/libxml2/include
        ${PROJECT_SOURCE_DIR}/lib/hdf5-1.8.12/include
)

link_directories(${PROJECT_SOURCE_DIR}/lib/gsl-1.14/lib)
link_directories(${PROJECT_SOURCE_DIR}/lib/libxml2/lib)
link_directories(${PROJECT_SOURCE_DIR}/lib/hdf5-1.8.12/lib)

target_link_libraries(sacess
        # ${LIBS}
        # ${CLIBS}
        # ${FLIBS}
)