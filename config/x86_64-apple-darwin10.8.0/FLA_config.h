/* build/FLA_config.h.  Generated from FLA_config.h.in by configure.  */
/* build/FLA_config.h.in.  Generated from configure.ac by autoheader.  */

/* Sets the default blocksize in the k dimension. */
/* #undef FLA_DEFAULT_K_BLOCKSIZE */

/* Sets the default blocksize in the m dimension. */
/* #undef FLA_DEFAULT_M_BLOCKSIZE */

/* Sets the default blocksize in the n dimension. */
/* #undef FLA_DEFAULT_N_BLOCKSIZE */

/* Determines whether to use control trees to select a reasonable FLAME
   variant and blocksize when level-3 BLAS front-ends are invoked. */
/* #undef FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES */

/* Determines whether to define bli_malloc() in terms of FLA_malloc(). */
#define FLA_ENABLE_BLIS_USE_OF_FLA_MALLOC 1

/* Determines whether a built-in implementation of the BLAS is compiled. */
/* #undef FLA_ENABLE_BUILTIN_BLAS */

/* Determines whether to enable CBLAS interfaces instead of Fortran-77
   interfaces to the BLAS. */
/* #undef FLA_ENABLE_CBLAS_INTERFACES */

/* Determines whether to enable external LAPACK for small subproblems. */
/* #undef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS */

/* Determines whether to enable interfaces to external LAPACK routines. */
/* #undef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES */

/* Determines whether to enable interfaces to internal/low-level libgoto
   symbols. */
/* #undef FLA_ENABLE_GOTO_INTERFACES */

/* Determines whether GPU-specific blocks of code should be compiled. */
/* #undef FLA_ENABLE_GPU */

/* Determines whether to enable internal runtime consistency checks of
   function parameters and return values. */
#define FLA_ENABLE_INTERNAL_ERROR_CHECKING 1

/* Determines whether the LAPACK compatibility layer is included in libflame.
   */
/* #undef FLA_ENABLE_LAPACK2FLAME */

/* Determines whether to enable code that will increase FLA_Obj leading
   dimensions to ensure that matrix columns adhere to the alignment specified
   by FLA_MEMORY_ALIGNMENT_BOUNDARY. */
/* #undef FLA_ENABLE_LDIM_ALIGNMENT */

/* Determines whether memory is aligned to user-requested boundaries. */
/* #undef FLA_ENABLE_MEMORY_ALIGNMENT */

/* Determines whether to enable the FLA_malloc()/FLA_free() memory counter by
   default. */
/* #undef FLA_ENABLE_MEMORY_LEAK_COUNTER */

/* Determines whether thread-specific blocks of code should be compiled. */
/* #undef FLA_ENABLE_MULTITHREADING */

/* Determines whether to enable various segments of code identified as
   providing non-critical functionality. */
#define FLA_ENABLE_NON_CRITICAL_CODE 1

/* Determines whether to define a portable FLA_Clock() in terms of
   clock_gettime() or gettimeofday() from time.h. */
#define FLA_ENABLE_PORTABLE_TIMER 1

/* Determines whether SCC-specific blocks of code should be compiled. */
/* #undef FLA_ENABLE_SCC */

/* Determines whether SuperMatrix-specific blocks of code should be compiled.
   */
/* #undef FLA_ENABLE_SUPERMATRIX */

/* Determines whether blocks of code specific to Texas Instruments' DSP. */
/* #undef FLA_ENABLE_TIDSP */

/* Determines whether vector intrinsics are used in certain low-level
   functions. */
/* #undef FLA_ENABLE_VECTOR_INTRINSICS */

/* Determines whether to modify various segments of code needed for
   integrating libflame into Windows. */
/* #undef FLA_ENABLE_WINDOWS_BUILD */

/* Encodes the default level of internal error checking chosen at
   configure-time. */
#define FLA_INTERNAL_ERROR_CHECKING_LEVEL 2

/* Sets the byte boundary used to align the starting address of all memory
   allocated dynamically through libflame. */
/* #undef FLA_MEMORY_ALIGNMENT_BOUNDARY */

/* Encodes the type of multithreading chosen at configure-time. */
#define FLA_MULTITHREADING_MODEL 0

/* Determines whether clock_gettime() was present on the system (via time.h).
   */
/* #undef FLA_PORTABLE_TIMER_IS_CLOCK_GETTIME */

/* Determines whether gettimeofday() was present on the system (via time.h).
   */
#define FLA_PORTABLE_TIMER_IS_GETTIMEOFDAY 1

/* Determines whether a timer was found at all. */
/* #undef FLA_PORTABLE_TIMER_IS_UNKNOWN */

/* Encodes the type of vector intrinsics requested at configure-time. */
#define FLA_VECTOR_INTRINSIC_TYPE 0

/* Define to 1 if you have the <assert.h> header file. */
#define HAVE_ASSERT_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the <ia64intrin.h> header file. */
/* #undef HAVE_IA64INTRIN_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <signal.h> header file. */
#define HAVE_SIGNAL_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "flame@cs.utexas.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "libflame"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "libflame r`cat revision`"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "libflame"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.cs.utexas.edu/users/flame/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "r`cat revision`"

/* Define to 1 if the C compiler supports function prototypes. */
#define PROTOTYPES 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Enables ANSI C, POSIX.1, POSIX.2, BSD, SVID, X/Open, and GNU extensions to
   the C language. */
#define _GNU_SOURCE 1

/* Define like PROTOTYPES; this can be used by system headers. */
#define __PROTOTYPES 1

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
