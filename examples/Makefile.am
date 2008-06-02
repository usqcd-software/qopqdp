ALL_CFLAGS = -I$(top_srcdir)/include -I$(top_builddir)/include -I$(QDPDIR)/include -I$(QLADIR)/include -I$(QIODIR)/include -I$(QMPDIR)/include
ALL_LDADD = -L$(top_builddir)/lib -L$(QDPDIR)/lib  -L$(QLADIR)/lib -L$(QIODIR)/lib -L$(QMPDIR)/lib -lqdp_common -lqdp_int -lqla_int -lqla_random -lqio -llime -lqmp -lm

F3_CFLAGS = -DQOP_Precision=1 -DQDP_Precision=1 -DQDP_Nc=3 $(ALL_CFLAGS)
F3_LDADD = -lqopqdp -lqdp_f -lqdp_f3 -lqdp_df3 -lqla_f -lqla_f3 -lqla_d -lqla_d3 -lqla_df -lqla_df3 $(ALL_LDADD)

D3_CFLAGS = -DQOP_Precision=2 -DQDP_Precision=2 -DQDP_Nc=3 $(ALL_CFLAGS)
D3_LDADD = -lqopqdp -lqdp_d -lqdp_d3 -lqdp_df3 -lqla_d -lqla_d3 -lqla_q -lqla_q3 -lqla_dq -lqla_dq3 $(ALL_LDADD)

check_PROGRAMS = \
  imp_gauge_force_bench_f3 imp_gauge_force_bench_d3 \
  asqtad_test_f3 asqtad_test_d3 \
  asqtad_invert_bench_f3 asqtad_invert_bench_d3 \
  asqtad_force_bench_f3 asqtad_force_bench_d3 \
  asqtad_eigcg_f3 asqtad_eigcg_d3 \
  hisq_invert_bench_f3 hisq_invert_bench_d3 \
  wilson_test_f3 wilson_test_d3 \
  wilson_invert_bench_f3 wilson_invert_bench_d3 \
  dw_invert_bench_f3 dw_invert_bench_d3

imp_gauge_force_bench_f3_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_f3_CFLAGS = $(F3_CFLAGS)
imp_gauge_force_bench_f3_LDADD = $(F3_LDADD)

imp_gauge_force_bench_d3_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_d3_CFLAGS = $(D3_CFLAGS)
imp_gauge_force_bench_d3_LDADD = $(D3_LDADD)

asqtad_test_f3_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_f3_CFLAGS = $(F3_CFLAGS)
asqtad_test_f3_LDADD = $(F3_LDADD)

asqtad_invert_bench_f3_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_f3_CFLAGS = $(F3_CFLAGS)
asqtad_invert_bench_f3_LDADD = $(F3_LDADD)

asqtad_force_bench_f3_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_f3_CFLAGS = $(F3_CFLAGS)
asqtad_force_bench_f3_LDADD = $(F3_LDADD)

asqtad_eigcg_f3_SOURCES = test_common.h test_common.c asqtad_eigcg.c
asqtad_eigcg_f3_CFLAGS = $(F3_CFLAGS)
asqtad_eigcg_f3_LDADD = $(F3_LDADD)

asqtad_test_d3_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_d3_CFLAGS = $(D3_CFLAGS)
asqtad_test_d3_LDADD = $(D3_LDADD)

asqtad_invert_bench_d3_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_d3_CFLAGS = $(D3_CFLAGS)
asqtad_invert_bench_d3_LDADD = $(D3_LDADD)

asqtad_force_bench_d3_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_d3_CFLAGS = $(D3_CFLAGS)
asqtad_force_bench_d3_LDADD = $(D3_LDADD)

asqtad_eigcg_d3_SOURCES = test_common.h test_common.c asqtad_eigcg.c
asqtad_eigcg_d3_CFLAGS = $(D3_CFLAGS)
asqtad_eigcg_d3_LDADD = $(D3_LDADD)


hisq_invert_bench_f3_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_f3_CFLAGS = $(F3_CFLAGS)
hisq_invert_bench_f3_LDADD = $(F3_LDADD)

hisq_invert_bench_d3_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_d3_CFLAGS = $(D3_CFLAGS)
hisq_invert_bench_d3_LDADD = $(D3_LDADD)


wilson_test_f3_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_f3_CFLAGS = $(F3_CFLAGS)
wilson_test_f3_LDADD = $(F3_LDADD)

wilson_invert_bench_f3_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_f3_CFLAGS = $(F3_CFLAGS)
wilson_invert_bench_f3_LDADD = $(F3_LDADD)

wilson_test_d3_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_d3_CFLAGS = $(D3_CFLAGS)
wilson_test_d3_LDADD = $(D3_LDADD)

wilson_invert_bench_d3_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_d3_CFLAGS = $(D3_CFLAGS)
wilson_invert_bench_d3_LDADD = $(D3_LDADD)


dw_invert_bench_f3_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_f3_CFLAGS = $(F3_CFLAGS)
dw_invert_bench_f3_LDADD = $(F3_LDADD)

dw_invert_bench_d3_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_d3_CFLAGS = $(D3_CFLAGS)
dw_invert_bench_d3_LDADD = $(D3_LDADD)
