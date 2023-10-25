# DO NOT CHANGE ANYTHING HERE!!! ##################################
# (unless you know *exactly* what you're doing...) 

#all objects
AOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o savestate.o fly_sa.o ../util/error.o ../util/ioTools.o solvers.o score.o # ../util/distributions.o ../util/random.o ../util/dSFMT.o ../util/dSFMT_str_state.o

#fly_sa objects
FOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
         ../util/error.o ../util/distributions.o ../util/random.o ../util/ioTools.o solvers.o score.o ../util/dSFMT.o ../util/dSFMT_str_state.o
# serial code
FSOBJ = fly_sa.o savestate.o # moves.o ../util/lsa.o
# parallel code
# FPOBJ = moves-mpi.o fly_sa-mpi.o ../util/lsa-mpi.o savestate-mpi.o


#printscore objects
POBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
       ../util/error.o ../util/distributions.o ../util/random.o ../util/ioTools.o solvers.o score.o printscore.o ../util/dSFMT.o ../util/dSFMT_str_state.o

#unfold objects
UOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
	 ../util/error.o ../util/distributions.o ../util/random.o ../util/ioTools.o solvers.o score.o unfold.o ../util/dSFMT.o ../util/dSFMT_str_state.o

#scramble objects
SOBJ = zygotic.o fly_io.o maternal.o integrate.o translate.o \
	 ../util/error.o ../util/distributions.o ../util/random.o ../util/ioTools.o solvers.o score.o scramble.o ../util/dSFMT.o ../util/dSFMT_str_state.o

SOURCES = `ls *.c`

#Below here are the rules for building things

all: $(FLYEXECS)

# special cases: dependencies and flags for individual .c files

fly_sa.o: fly_sa.c
	$(CC) -c $(CFLAGS) $(VFLAGS) fly_sa.c

printscore.o: printscore.c
	$(CC) -c $(CFLAGS) $(VFLAGS) printscore.c

scramble.o: scramble.c
	$(CC) -c $(CFLAGS) $(VFLAGS) scramble.c

unfold.o: unfold.c
	$(CC) -c $(CFLAGS) $(VFLAGS) unfold.c

zygotic.o: zygotic.c
	$(CC) -c $(CFLAGS) zygotic.c

objects: $(AOBJ)	


# parallel stuff

#fly_sa-mpi.o: fly_sa.c
#	$(MPICC) -c -o fly_sa-mpi.o $(MPIFLAGS) $(CFLAGS) $(VFLAGS) fly_sa.c

#lsa-mpi.o: ../util/lsa.c
#	$(MPICC) -c -o ../util/lsa-mpi.o $(MPIFLAGS) $(CFLAGS) ../util/lsa.c

#moves-mpi.o: moves.c 
#	$(MPICC) -c -o moves-mpi.o $(MPIFLAGS) $(CFLAGS) moves.c

#savestate-mpi.o: savestate.c
#	$(MPICC) -c -o savestate-mpi.o $(MPIFLAGS) $(CFLAGS) savestate.c

# executable targets: serial ...

fly_sa: $(FOBJ) $(FSOBJ)
	$(CC) -o fly_sa $(CFLAGS) $(LDFLAGS) $(FOBJ) $(FSOBJ) $(FLIBS) 

printscore: $(POBJ)
	$(CC) -o printscore $(CFLAGS) $(LDFLAGS) $(POBJ) $(LIBS) 

unfold: $(UOBJ)
	$(CC) -o unfold $(CFLAGS) $(LDFLAGS) $(UOBJ) $(LIBS) 

scramble: $(SOBJ)
	$(CC) -o scramble $(CFLAGS) $(LDFLAGS) $(SOBJ) $(LIBS) 

# ... and parallel

#fly_sa.mpi: $(FOBJ) $(FPOBJ)
#	$(MPICC) -o fly_sa.mpi $(CFLAGS) $(LDFLAGS) $(FOBJ) $(FPOBJ) $(FLIBS)

# ... and here are the cleanup and make deps rules

clean:
	rm -f *.o core*

Makefile: ${FRC}
	rm -f $@
	cp basic.mk $@
	echo "#Automatically generated dependencies list#" >> $@
	${CC} $(INCLUDES) -M ${SOURCES} >> $@
	chmod -w $@

#Automatically generated dependencies list#
fly_io.o: fly_io.c /usr/include/stdc-predef.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/mathLib.h \
 fly_io.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/ctype.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h maternal.h
fly_sa.o: fly_sa.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/limits.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/syslimits.h \
 /usr/include/limits.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/include/x86_64-linux-gnu/bits/posix1_lim.h \
 /usr/include/x86_64-linux-gnu/bits/local_lim.h \
 /usr/include/linux/limits.h \
 /usr/include/x86_64-linux-gnu/bits/pthread_stack_min-dynamic.h \
 /usr/include/x86_64-linux-gnu/bits/pthread_stack_min.h \
 /usr/include/x86_64-linux-gnu/bits/posix2_lim.h /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/include/unistd.h \
 /usr/include/x86_64-linux-gnu/bits/posix_opt.h \
 /usr/include/x86_64-linux-gnu/bits/environments.h \
 /usr/include/x86_64-linux-gnu/bits/confname.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_posix.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_core.h \
 /usr/include/x86_64-linux-gnu/bits/unistd_ext.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 maternal.h zygotic.h /usr/include/time.h \
 /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h score.h \
 fly_io.h /usr/include/ctype.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/math.h /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h \
 /usr/include/x86_64-linux-gnu/sys/time.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/sa.h \
 fly_sa.h
integrate.o: integrate.c /usr/include/stdc-predef.h /usr/include/stdio.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/include/x86_64-linux-gnu/sys/times.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 maternal.h zygotic.h /usr/include/time.h \
 /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h score.h \
 fly_io.h /usr/include/ctype.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/x86_64-linux-gnu/sys/time.h
maternal.o: maternal.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h /usr/include/ctype.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/string.h /usr/include/strings.h maternal.h \
 /usr/include/stdio.h /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 zygotic.h /usr/include/time.h /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h score.h \
 fly_io.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/x86_64-linux-gnu/sys/time.h
printscore.o: printscore.c /usr/include/stdc-predef.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/include/unistd.h \
 /usr/include/x86_64-linux-gnu/bits/posix_opt.h \
 /usr/include/x86_64-linux-gnu/bits/environments.h \
 /usr/include/x86_64-linux-gnu/bits/confname.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_posix.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_core.h \
 /usr/include/x86_64-linux-gnu/bits/unistd_ext.h integrate.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 maternal.h zygotic.h /usr/include/time.h \
 /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h score.h \
 fly_io.h /usr/include/ctype.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/x86_64-linux-gnu/sys/time.h maternal.h score.h solvers.h \
 zygotic.h
savestate.o: savestate.c /usr/include/stdc-predef.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/unistd.h \
 /usr/include/x86_64-linux-gnu/bits/posix_opt.h \
 /usr/include/x86_64-linux-gnu/bits/environments.h \
 /usr/include/x86_64-linux-gnu/bits/confname.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_posix.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_core.h \
 /usr/include/x86_64-linux-gnu/bits/unistd_ext.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/random.h \
 fly_io.h /usr/include/ctype.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h maternal.h
score.o: score.c /usr/include/stdc-predef.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/ctype.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/string.h /usr/include/strings.h score.h maternal.h fly_io.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 zygotic.h /usr/include/time.h /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/x86_64-linux-gnu/sys/time.h
scramble.o: scramble.c /usr/include/stdc-predef.h /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/include/unistd.h \
 /usr/include/x86_64-linux-gnu/bits/posix_opt.h \
 /usr/include/x86_64-linux-gnu/bits/environments.h \
 /usr/include/x86_64-linux-gnu/bits/confname.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_posix.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_core.h \
 /usr/include/x86_64-linux-gnu/bits/unistd_ext.h maternal.h score.h \
 maternal.h fly_io.h /usr/include/ctype.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 zygotic.h /usr/include/time.h /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h score.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/math.h /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h \
 /usr/include/x86_64-linux-gnu/sys/time.h zygotic.h ../util/random.h
solvers.o: solvers.c /usr/include/stdc-predef.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h /usr/include/stdlib.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/lib/gcc/x86_64-linux-gnu/11/include/limits.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/syslimits.h \
 /usr/include/limits.h /usr/include/x86_64-linux-gnu/bits/posix1_lim.h \
 /usr/include/x86_64-linux-gnu/bits/local_lim.h \
 /usr/include/linux/limits.h \
 /usr/include/x86_64-linux-gnu/bits/pthread_stack_min-dynamic.h \
 /usr/include/x86_64-linux-gnu/bits/pthread_stack_min.h \
 /usr/include/x86_64-linux-gnu/bits/posix2_lim.h solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 maternal.h maternal.h zygotic.h /usr/include/time.h \
 /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 solvers.h zygotic.h score.h fly_io.h /usr/include/ctype.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/x86_64-linux-gnu/sys/time.h
translate.o: translate.c /usr/include/stdc-predef.h /usr/include/ctype.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/string.h \
 /usr/include/strings.h translate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h maternal.h \
 /usr/include/stdio.h /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h score.h zygotic.h \
 /usr/include/time.h /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/math.h /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h \
 /usr/include/x86_64-linux-gnu/sys/time.h fly_io.h
unfold.o: unfold.c /usr/include/stdc-predef.h /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/stdio.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h /usr/include/string.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/strings.h /usr/include/unistd.h \
 /usr/include/x86_64-linux-gnu/bits/posix_opt.h \
 /usr/include/x86_64-linux-gnu/bits/environments.h \
 /usr/include/x86_64-linux-gnu/bits/confname.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_posix.h \
 /usr/include/x86_64-linux-gnu/bits/getopt_core.h \
 /usr/include/x86_64-linux-gnu/bits/unistd_ext.h /usr/include/time.h \
 /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h integrate.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 maternal.h zygotic.h score.h fly_io.h /usr/include/ctype.h integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/math.h /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h \
 /usr/include/x86_64-linux-gnu/sys/time.h maternal.h solvers.h zygotic.h \
 score.h fly_io.h
zygotic.o: zygotic.c /usr/include/stdc-predef.h /usr/include/ctype.h \
 /usr/include/features.h /usr/include/features-time64.h \
 /usr/include/x86_64-linux-gnu/bits/wordsize.h \
 /usr/include/x86_64-linux-gnu/bits/timesize.h \
 /usr/include/x86_64-linux-gnu/sys/cdefs.h \
 /usr/include/x86_64-linux-gnu/bits/long-double.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs.h \
 /usr/include/x86_64-linux-gnu/gnu/stubs-64.h \
 /usr/include/x86_64-linux-gnu/bits/types.h \
 /usr/include/x86_64-linux-gnu/bits/typesizes.h \
 /usr/include/x86_64-linux-gnu/bits/time64.h \
 /usr/include/x86_64-linux-gnu/bits/endian.h \
 /usr/include/x86_64-linux-gnu/bits/endianness.h \
 /usr/include/x86_64-linux-gnu/bits/types/locale_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__locale_t.h \
 /usr/include/stdio.h \
 /usr/include/x86_64-linux-gnu/bits/libc-header-start.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stddef.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/stdarg.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__mbstate_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__fpos64_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/FILE.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_FILE.h \
 /usr/include/x86_64-linux-gnu/bits/stdio_lim.h \
 /usr/include/x86_64-linux-gnu/bits/floatn.h \
 /usr/include/x86_64-linux-gnu/bits/floatn-common.h /usr/include/stdlib.h \
 /usr/include/x86_64-linux-gnu/bits/waitflags.h \
 /usr/include/x86_64-linux-gnu/bits/waitstatus.h \
 /usr/include/x86_64-linux-gnu/sys/types.h \
 /usr/include/x86_64-linux-gnu/bits/types/clock_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/clockid_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/time_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/timer_t.h \
 /usr/include/x86_64-linux-gnu/bits/stdint-intn.h /usr/include/endian.h \
 /usr/include/x86_64-linux-gnu/bits/byteswap.h \
 /usr/include/x86_64-linux-gnu/bits/uintn-identity.h \
 /usr/include/x86_64-linux-gnu/sys/select.h \
 /usr/include/x86_64-linux-gnu/bits/select.h \
 /usr/include/x86_64-linux-gnu/bits/types/sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/__sigset_t.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timeval.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_timespec.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes.h \
 /usr/include/x86_64-linux-gnu/bits/thread-shared-types.h \
 /usr/include/x86_64-linux-gnu/bits/pthreadtypes-arch.h \
 /usr/include/x86_64-linux-gnu/bits/atomic_wide_counter.h \
 /usr/include/x86_64-linux-gnu/bits/struct_mutex.h \
 /usr/include/x86_64-linux-gnu/bits/struct_rwlock.h /usr/include/alloca.h \
 /usr/include/x86_64-linux-gnu/bits/stdlib-float.h /usr/include/string.h \
 /usr/include/strings.h /usr/include/math.h \
 /usr/include/x86_64-linux-gnu/bits/math-vector.h \
 /usr/include/x86_64-linux-gnu/bits/libm-simd-decl-stubs.h \
 /usr/include/x86_64-linux-gnu/bits/flt-eval-method.h \
 /usr/include/x86_64-linux-gnu/bits/fp-logb.h \
 /usr/include/x86_64-linux-gnu/bits/fp-fast.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls-helper-functions.h \
 /usr/include/x86_64-linux-gnu/bits/mathcalls.h maternal.h solvers.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_nvector.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_config.h \
 /usr/lib/gcc/x86_64-linux-gnu/11/include/float.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_direct.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_band.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/cvode/cvode_spils.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_iterative.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_spgmr.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/nvector/nvector_serial.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_dense.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/sundials-2.5.0/include/sundials/sundials_math.h \
 integrate.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_errno.h \
 /usr/include/errno.h /usr/include/x86_64-linux-gnu/bits/errno.h \
 /usr/include/linux/errno.h /usr/include/x86_64-linux-gnu/asm/errno.h \
 /usr/include/asm-generic/errno.h /usr/include/asm-generic/errno-base.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_types.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_spline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_interp.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/lib64/gsl-1.14/include/gsl/gsl_inline.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/error.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/global.h \
 zygotic.h /usr/include/time.h /usr/include/x86_64-linux-gnu/bits/time.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_tm.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_itimerspec.h \
 /usr/include/x86_64-linux-gnu/sys/resource.h \
 /usr/include/x86_64-linux-gnu/bits/resource.h \
 /usr/include/x86_64-linux-gnu/bits/types/struct_rusage.h score.h \
 fly_io.h \
 /home/david/biopredynlib/source/benchmarks/systemsBiology/BioPredyn/B6/util/ioTools.h \
 /usr/include/x86_64-linux-gnu/sys/time.h
