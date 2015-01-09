AC_INIT([mapmap], [0.13.0], [astro@joergdietrich.com])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE(foreign)

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET


AC_SUBST(CFLAGS,$CFLAGS)
AC_SUBST(LDFLAGS, $LDFLAGS)

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h libintl.h limits.h malloc.h stddef.h stdlib.h string.h unistd.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday memset strcasecmp strdup])

AC_CHECK_PROG(HAVE_GSL, gsl-config, [yes], [no])
CFLAGS="$CFLAGS `gsl-config --cflags`"
LDFLAGS="$LDFLAGS `gsl-config --libs`"
AC_CHECK_LIB(m,cos)
AC_CHECK_LIB(gslcblas,cblas_dgemm, [], \
                 [echo ; echo Linking against gsl failed.;
                 echo Make sure that GSL with CBLAS is installed.; echo ;
                 exit 1])
AC_CHECK_LIB(gsl,gsl_blas_dgemm, [], \
                 [echo ; echo Linking against gsl failed.;
                 echo Make sure that GSL with CBLAS is installed.; echo ;
                 exit 1])
AC_CHECK_HEADER(fitsio.h, [], \
			  [echo ; echo fitsio.h not found.;\
			  echo Make sure that cfitsio is installed.; echo ;
			  exit 1], \
			  [])
AC_CHECK_LIB(cfitsio, ffopen, [], \
		      [echo ; echo Linking against cfitsio failed.;
		      echo Make sure that cfitsio is installed.; echo ;
		      exit 1], \
		      [-lm])

AC_MSG_CHECKING(whether to enable parallel algorithms with OpenMP)
AC_ARG_ENABLE([openmp], [AC_HELP_STRING([--enable-openmp],
			[enable OpenMP [default=no] ])],
			if test $enableval = yes; then
			   [AC_MSG_RESULT(yes)]
			   [AC_DEFINE([USING_OPENMP], 1, [using openmp])]
			   USING_OPENMP=yes
			   [CFLAGS="$CFLAGS -fopenmp "]
			   [LDFLAGS="$LDFLAGS -lgomp "]
			fi
			,
			[AC_MSG_RESULT(no)]
			USING_OPENMP=no
			[LDFLAGS="$LDFLAGS"]
			)

AC_SUBST(OPENMP_LDFLAGS)

AC_MSG_CHECKING(whether to enable Portable Linux Processor Affinity)
AC_ARG_ENABLE([plpa], [AC_HELP_STRING([--enable-plpa],
		      [enable PLPA [default=no] ])],
		      if test $enableval = yes; then
			 [AC_MSG_RESULT(yes)]
			 [AC_DEFINE([USING_PLPA], 1, [using plpa])]
			 USING_PLPA=yes
			 [PLPA_INCLUDED(plpa-1.1)]
			 [PLPA_SET_SYMBOL_PREFIX(mplpa_)]
			 [PLPA_INIT(plpa_happy=yes, plpa_happy=no)]
			 if test $plpa_happy = no; then
			    echo PLPA is not happy
			    exit 1
			 fi
		      fi
		      ,
		      [AC_MSG_RESULT(no)]
		      [PLPA_DO_AM_CONDITIONALS]
		      USING_PLPA=no
		      ) 

AM_CONDITIONAL([PLPA_OPT], [test "$plpa_happy" = yes])
AC_CONFIG_FILES(Makefile \
          lib/Makefile \
          src/Makefile)
AC_OUTPUT


