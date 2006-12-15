AM_CFLAGS = -I@top_srcdir@/include

lib_LIBRARIES = libqopqdp.a

libqopqdp_a_SOURCES = \
  common.c common_f.c common_d.c \
  symanzik_1loop_gauge_force_f.c symanzik_1loop_gauge_force_d.c \
  asqtad_invert.c asqtad_invert_f.c asqtad_invert_d.c \
  asqtad_force.c asqtad_force_f.c asqtad_force_d.c \
  asqtad_force_fnmat_f.c asqtad_force_fnmat_d.c \
  wilson_invert.c wilson_invert_f.c wilson_invert_d.c \
  dw_invert.c dw_invert_f.c dw_invert_d.c \
  invert_cg_V_f.c invert_cg_V_d.c \
  invert_cg_D_f.c invert_cg_D_d.c \
  invert_cg_vD_f.c invert_cg_vD_d.c

common_f.o: common_p.c
common_d.o: common_p.c
symanzik_1loop_gauge_force_f.c: symanzik_1loop_gauge_force_p.c
symanzik_1loop_gauge_force_d.c: symanzik_1loop_gauge_force_p.c
asqtad_invert_f.o: asqtad_invert_p.c
asqtad_invert_d.o: asqtad_invert_p.c
asqtad_force_f.o: asqtad_force_p.c
asqtad_force_d.o: asqtad_force_p.c
asqtad_force_fnmat_f.o: asqtad_force_fnmat_p.c
asqtad_force_fnmat_d.o: asqtad_force_fnmat_p.c
wilson_invert_f.o: wilson_invert_p.c
wilson_invert_d.o: wilson_invert_p.c
dw_invert_f.o: dw_invert_p.c
dw_invert_d.o: dw_invert_p.c
invert_cg_V_f.o: invert_cg_p.c
invert_cg_V_d.o: invert_cg_p.c
invert_cg_D_f.o: invert_cg_p.c
invert_cg_D_d.o: invert_cg_p.c
invert_cg_vD_f.o: invert_cg_p.c
invert_cg_vD_d.o: invert_cg_p.c

EXTRA_DIST = \
  asqtad_action.h \
  common_p.c \
  symanzik_1loop_gauge_force_p.c \
  asqtad_invert_p.c \
  asqtad_force_p.c \
  asqtad_force_fnmat_p.c \
  wilson_invert_p.c \
  dw_invert_p.c \
  invert_cg_p.c



#lib_LIBRARIES = libqopqdp_common.a libqopqdp_f.a libqopqdp_d.a

#LIBCOMMONSOURCE = \
#  ../include/qop.h \
#  ../include/qop_qdp.h \
#  ../include/qop_internal.h \
#  common.c \
#  asqtad.c \
#  wilson.c

#LIBPSOURCE = \
#  ../include/qop.h \
#  ../include/qop_qdp.h \
#  ../include/qop_internal.h \
#  asqtadp.c \
#  wilsonp.c

#libqopqdp_common_a_SOURCES = $(LIBCOMMONSOURCE)
#libqopqdp_common_a_CFLAGS = $(AM_CFLAGS)

#libqopqdp_f_a_SOURCES = $(LIBPSOURCE)
#libqopqdp_f_a_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=1 -DQDP_Precision=1

#libqopqdp_d_a_SOURCES = $(LIBPSOURCE)
#libqopqdp_d_a_CFLAGS = $(AM_CFLAGS) -DQOP_Precision=2 -DQDP_Precision=2




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
