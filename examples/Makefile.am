ALL_CFLAGS = -I$(top_srcdir)/include -I$(top_builddir)/include -I$(QDPDIR)/include -I$(QLADIR)/include -I$(QIODIR)/include -I$(QMPDIR)/include
ALL_LDADD = -L$(top_builddir)/lib -L$(QDPDIR)/lib  -L$(QLADIR)/lib -L$(QIODIR)/lib -L$(QMPDIR)/lib -lqopqdp_common -lqdp_common -lqdp_int -lqla_int -lqla_random -lqio -llime -lqmp -lm

F3_CFLAGS = -DQOP_Precision=1 -DQDP_Precision=1 -DQDP_Nc=3 $(ALL_CFLAGS)
F3_LDADD = -lqopqdp_f -lqdp_f -lqdp_f3 -lqla_f -lqla_f3 -lqla_d -lqla_d3 -lqla_df -lqla_df3 $(ALL_LDADD)

check_PROGRAMS = asqtad_test_f3

asqtad_test_f3_SOURCES = asqtad_test.c
asqtad_test_f3_CFLAGS = $(F3_CFLAGS)
asqtad_test_f3_LDADD = $(F3_LDADD)

