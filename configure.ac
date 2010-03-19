AC_INIT(qopqdp, 0.16.2, [detar at physics.utah.edu])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR(include/qop.h)
AC_CONFIG_SRCDIR(lib/common.c)
AM_CONFIG_HEADER(include/qop_config.h)

dnl --enable-blas
AC_ARG_ENABLE(
  blas,
  AC_HELP_STRING([--enable-blas],
    [enable routines that use blas (set LDFLAGS and LIBS to link to blas, see --enable-underscores)]),
  if test "X${enableval}X" == "XyesX"; then
    AC_DEFINE(HAVE_BLAS, [], [have blas])
  fi
)

dnl --enable-lapack
AC_ARG_ENABLE(
  lapack,
  AC_HELP_STRING([--enable-lapack],
    [enable routines that use lapack and blas (set LDFLAGS and LIBS to link to lapack and blas, see --enable-underscores)]),
  if test "X${enableval}X" == "XyesX"; then
    AC_DEFINE(HAVE_BLAS, [], [have blas])
    AC_DEFINE(HAVE_LAPACK, [], [have lapack])
  fi
)

dnl --enable-underscores
AC_ARG_ENABLE(
  underscores,
  AC_HELP_STRING([--enable-underscores],
    [append underscores to lapack and blas routines]),
  if test "X${enableval}X" == "XyesX"; then
    AC_DEFINE(USE_UNDERSCORES, [], [use underscores])
  fi
)

AC_ARG_WITH(
  qdp,
  AC_HELP_STRING([--with-qdp=<QDP dir>],
    [location of QDP libraries]),
  [QDPDIR="$with_qdp"; WITH_QDP=1]
)
AC_SUBST(QDPDIR)

AC_ARG_WITH(
  qla,
  AC_HELP_STRING([--with-qla=<QLA dir>],
    [location of QLA libraries]),
  [QLADIR="$with_qla"; WITH_QLA=1]
)
AC_SUBST(QLADIR)

AC_ARG_WITH(
  qio,
  AC_HELP_STRING([--with-qio=<QIO dir>],
    [location of QIO libraries]),
  [QIODIR="$with_qio"; WITH_QIO=1]
)
AC_SUBST(QIODIR)

AC_ARG_WITH(
  qmp,
  AC_HELP_STRING([--with-qmp=<QMP dir>],
    [location of QMP libraries]),
  [QMPDIR="$with_qmp"; WITH_QMP=1]
)
AC_SUBST(QMPDIR)

AC_SUBST(QDP_CFLAGS, "-I${QDPDIR}/include")
AC_SUBST(QDP_LDFLAGS, "-L${QDPDIR}/lib")
#AC_SUBST(QDP_LIBS, "-lqdp_fn -lqdp_f -lqdp_int -lqdp_common")
AC_SUBST(QDP_LIBS, "-lqdp_f3 -lqdp_f -lqdp_d3 -lqdp_d -lqdp_int -lqdp_common")

AC_SUBST(QLA_CFLAGS, "-I${QLADIR}/include")
AC_SUBST(QLA_LDFLAGS, "-L${QLADIR}/lib")
#AC_SUBST(QLA_LIBS, "-lqla_dfn -lqla_df -lqla_dn -lqla_d -lqla_fn -lqla_f -lqla_int -lqla_random -lqla_cmath -lm")
AC_SUBST(QLA_LIBS, "-lqla_df3 -lqla_df -lqla_f3 -lqla_f -lqla_d3 -lqla_d -lqla_int -lqla_random -lqla_cmath -lm")

AC_SUBST(QIO_CFLAGS, "-I${QIODIR}/include")
AC_SUBST(QIO_LDFLAGS, "-L${QIODIR}/lib")
AC_SUBST(QIO_LIBS, "-lqio -llime")

# get QMP configuration
if test "X${QMPDIR}X" = "XX" ; then
  AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${PATH}])
else
  AC_PATH_PROG([QMP_CONFIG], [qmp-config], [], [${QMPDIR}/bin ${PATH}])
fi

if test "X${QMP_CONFIG}X" != "XX" ; then
  AC_MSG_NOTICE([Found QMP configuration program ${QMP_CONFIG}])
  # use the QMP CC unless it is explicitly set by the user
  if ! set | grep -q '^CC=' ; then
    AC_SUBST(CC,  "`${QMP_CONFIG} --cc`")
    AC_MSG_NOTICE([QMP compiler: ${CC}])
  fi
  AC_SUBST(QMP_CFLAGS,  "`${QMP_CONFIG} --cflags`")
  AC_MSG_NOTICE([QMP compile flags: ${QMP_CFLAGS}])
  AC_SUBST(QMP_LDFLAGS, "`${QMP_CONFIG} --ldflags`")
  AC_MSG_NOTICE([QMP linking flags: ${QMP_LDFLAGS}])
  AC_SUBST(QMP_LIBS,    "`${QMP_CONFIG} --libs`")
  AC_MSG_NOTICE([QMP libraries flags: ${QMP_LIBS}])
else
  AC_MSG_WARN([QMP configuration program qmp-config not found.])
  AC_MSG_WARN([Set environment variables QMP_CFLAGS QMP_LDFAGS QMP_LIBS
    before configure])
fi

CFLAGS="$CFLAGS $QDP_CFLAGS $QLA_CFLAGS $QIO_CFLAGS $QMP_CFLAGS"
LDFLAGS="$LDFLAGS $QDP_LDFLAGS $QLA_LDFLAGS $QIO_LDFLAGS $QMP_LDFLAGS"
LIBS="$LIBS $QDP_LIBS $QLA_LIBS $QIO_LIBS $QMP_LIBS"

# Checks for programs.
AC_PROG_RANLIB
#AC_DISABLE_SHARED
#AC_PROG_LIBTOOL
AC_PROG_CC([c99 cc gcc])
AC_CHECK_TOOL(AR, ar, [ar])
AM_PROG_CC_C_O
#AC_PROG_CC_C99

# pass configuration to distcheck target
DISTCHECK_CONFIGURE_FLAGS="CC=\"$CC\" CFLAGS=\"$CFLAGS\" LDFLAGS=\"$LDFLAGS\" LIBS=\"$LIBS\""
if test -n "$WITH_QMP"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qmp=$QMPDIR"
fi
if test -n "$WITH_QIO"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qio=$QIODIR"
fi
if test -n "$WITH_QLA"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qla=$QLADIR"
fi
if test -n "$WITH_QDP"; then
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-qdp=$QDPDIR"
fi
AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(examples/Makefile)
AC_OUTPUT
