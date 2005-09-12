AM_CFLAGS = -I@top_srcdir@/include

lib_LIBRARIES = libqopqdp_common.a libqopqdp_f.a libqopqdp_d.a

LIBCOMMONSOURCE = \
  ../include/qop.h \
  ../include/qop_internal.h \
  common.c \
  asqtad.c \
  wilson.c

LIBPSOURCE = \
  ../include/qop.h \
  ../include/qop_internal.h \
  asqtadp.c \
  wilsonp.c

libqopqdp_common_a_SOURCES = $(LIBCOMMONSOURCE)
libqopqdp_common_a_CFLAGS = $(AM_CFLAGS)

libqopqdp_f_a_SOURCES = $(LIBPSOURCE)
libqopqdp_f_a_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=1 -DQDP_Precision=1

libqopqdp_d_a_SOURCES = $(LIBPSOURCE)
libqopqdp_d_a_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=2 -DQDP_Precision=2




#lib_LTLIBRARIES = libqopqdp.la
#noinst_LTLIBRARIES = libqopqdpf.la libqopqdpd.la

#libqopqdp_la_SOURCES = ../include/qop.h asqtad.c
#libqopqdp_la_LIBADD = libqopqdpf.la
#libqopqdp_la_LIBADD = libqopqdpf.la libqopqdpd.la 

#libqopqdpf_la_SOURCES = asqtadp.c
#libqopqdpf_la_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=1 -DQDP_Precision=1

#libqopqdpd_la_SOURCES = asqtadp.c
#libqopqdpd_la_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=2 -DQDP_Precision=2

#asqtadf.o: asqtad.c
#asqtadd.o: asqtad.c

#asqtadf_o_SOURCES = asqtad.c
#asqtadf_o_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=1 -DQDP_Precision=1

#asqtadd_o_SOURCES = asqtad.c
#asqtadd_o_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=2 -DQDP_Precision=2
