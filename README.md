(1) Compile library:

$ cd source;
$ make very clean;
$ make;
$ cd ..;

(2) Link the library and their dependencies. Example:

$ gcc -o main main_file.c source/bin/libAMIGObench.a -I./source/include/ -I./source/lib/libAMIGO/include/include_amigo -I./source/lib/libAMIGO/include/include_cvodes -lm -lhdf5 -O3 -cpp -DGNU -fPIC -no-pie  -L./source/lib/hdf5-1.8.12/lib  -lhdf5 -L./source/lib/BLAS -lblas -Dnullptr=0  -L./source/lib/gsl-1.14/lib -lgsl -fPIC -DEXPORT -lstdc++ -lpthread -lrt -lgfortran  -cpp -MMD -lm -ldl -lz
$ ./main

LOAD CIRCADIAN PROBLEM
cost CIRCADIAN  0.00000005523204201628
LOAD SSP PROBLEM
cost SSP  0.00000016927476889735
LOAD B2 PROBLEM
cost B2  383.42589345361125197087
