ALL_CFLAGS = -I$(top_srcdir)/include -I$(top_builddir)/include $(CFLAGS)
F2_CFLAGS = -DQOP_PrecisionInt=1 -DQOP_Nc=2 -DQDP_Precision=\'F\' -DQDP_Nc=2 $(ALL_CFLAGS)
D2_CFLAGS = -DQOP_PrecisionInt=2 -DQOP_Nc=2 -DQDP_Precision=\'D\' -DQDP_Nc=2 $(ALL_CFLAGS)
F3_CFLAGS = -DQOP_PrecisionInt=1 -DQOP_Nc=3 -DQDP_Precision=\'F\' -DQDP_Nc=3 $(ALL_CFLAGS)
D3_CFLAGS = -DQOP_PrecisionInt=2 -DQOP_Nc=3 -DQDP_Precision=\'D\' -DQDP_Nc=3 $(ALL_CFLAGS)
FN_CFLAGS = -DQOP_PrecisionInt=1 -DQOP_Nc=4 -DQOP_Colors=\'N\' -DQDP_Precision=\'F\' -DQDP_Nc=4 -DQDP_Colors=\'N\' $(ALL_CFLAGS)
DN_CFLAGS = -DQOP_PrecisionInt=2 -DQOP_Nc=4 -DQOP_Colors=\'N\' -DQDP_Precision=\'D\' -DQDP_Nc=4 -DQDP_Colors=\'N\' $(ALL_CFLAGS)
LIBS = -L$(top_builddir)/lib -lqopqdp @LIBS@

check_PROGRAMS =
if NCN
check_PROGRAMS += \
  imp_gauge_force_bench_fn imp_gauge_force_bench_dn \
  asqtad_test_fn asqtad_test_dn \
  asqtad_invert_bench_fn asqtad_invert_bench_dn \
  asqtad_force_bench_fn asqtad_force_bench_dn \
  hisq_invert_bench_fn hisq_invert_bench_dn \
  hisq_force_bench_fn hisq_force_bench_dn \
  wilson_test_fn wilson_test_dn \
  wilson_invert_bench_fn wilson_invert_bench_dn \
  wilson_force_bench_fn wilson_force_bench_dn \
  dw_invert_bench_fn dw_invert_bench_dn
endif
if NC2
check_PROGRAMS += \
  imp_gauge_force_bench_f2 imp_gauge_force_bench_d2 \
  asqtad_test_f2 asqtad_test_d2 \
  asqtad_invert_bench_f2 asqtad_invert_bench_d2 \
  asqtad_force_bench_f2 asqtad_force_bench_d2 \
  hisq_invert_bench_f2 hisq_invert_bench_d2 \
  hisq_force_bench_f2 hisq_force_bench_d2 \
  wilson_test_f2 wilson_test_d2 \
  wilson_invert_bench_f2 wilson_invert_bench_d2 \
  wilson_force_bench_f2 wilson_force_bench_d2 \
  dw_invert_bench_f2 dw_invert_bench_d2
endif
if NC3
check_PROGRAMS += \
  imp_gauge_force_bench_f3 imp_gauge_force_bench_d3 \
  asqtad_test_f3 asqtad_test_d3 \
  asqtad_invert_bench_f3 asqtad_invert_bench_d3 \
  asqtad_force_bench_f3 asqtad_force_bench_d3 \
  asqtad_eigcg_f3 asqtad_eigcg_d3 \
  hisq_invert_bench_f3 hisq_invert_bench_d3 \
  hisq_force_bench_f3 hisq_force_bench_d3 \
  wilson_test_f3 wilson_test_d3 \
  wilson_invert_bench_f3 wilson_invert_bench_d3 \
  wilson_force_bench_f3 wilson_force_bench_d3 \
  dw_invert_bench_f3 dw_invert_bench_d3 \
  test_reunit_utils_f3 test_reunit_utils_d3
endif

imp_gauge_force_bench_f2_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_f2_CFLAGS = $(F2_CFLAGS)

imp_gauge_force_bench_d2_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_d2_CFLAGS = $(D2_CFLAGS)

imp_gauge_force_bench_f3_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_f3_CFLAGS = $(F3_CFLAGS)

imp_gauge_force_bench_d3_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_d3_CFLAGS = $(D3_CFLAGS)

imp_gauge_force_bench_fn_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_fn_CFLAGS = $(FN_CFLAGS)

imp_gauge_force_bench_dn_SOURCES = test_common.h test_common.c imp_gauge_force_bench.c
imp_gauge_force_bench_dn_CFLAGS = $(DN_CFLAGS)


asqtad_test_f2_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_f2_CFLAGS = $(F2_CFLAGS)

asqtad_test_d2_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_d2_CFLAGS = $(D2_CFLAGS)

asqtad_test_f3_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_f3_CFLAGS = $(F3_CFLAGS)

asqtad_test_d3_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_d3_CFLAGS = $(D3_CFLAGS)

asqtad_test_fn_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_fn_CFLAGS = $(FN_CFLAGS)

asqtad_test_dn_SOURCES = test_common.h test_common.c asqtad_test.c
asqtad_test_dn_CFLAGS = $(DN_CFLAGS)


asqtad_invert_bench_f2_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_f2_CFLAGS = $(F2_CFLAGS)

asqtad_invert_bench_d2_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_d2_CFLAGS = $(D2_CFLAGS)

asqtad_invert_bench_f3_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_f3_CFLAGS = $(F3_CFLAGS)

asqtad_invert_bench_d3_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_d3_CFLAGS = $(D3_CFLAGS)

asqtad_invert_bench_fn_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_fn_CFLAGS = $(FN_CFLAGS)

asqtad_invert_bench_dn_SOURCES = test_common.h test_common.c asqtad_invert_bench.c
asqtad_invert_bench_dn_CFLAGS = $(DN_CFLAGS)


asqtad_force_bench_f2_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_f2_CFLAGS = $(F2_CFLAGS)

asqtad_force_bench_d2_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_d2_CFLAGS = $(D2_CFLAGS)

asqtad_force_bench_f3_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_f3_CFLAGS = $(F3_CFLAGS)

asqtad_force_bench_d3_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_d3_CFLAGS = $(D3_CFLAGS)

asqtad_force_bench_fn_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_fn_CFLAGS = $(FN_CFLAGS)

asqtad_force_bench_dn_SOURCES = test_common.h test_common.c asqtad_force_bench.c
asqtad_force_bench_dn_CFLAGS = $(DN_CFLAGS)


asqtad_eigcg_f3_SOURCES = test_common.h test_common.c asqtad_eigcg.c
asqtad_eigcg_f3_CFLAGS = $(F3_CFLAGS)

asqtad_eigcg_d3_SOURCES = test_common.h test_common.c asqtad_eigcg.c
asqtad_eigcg_d3_CFLAGS = $(D3_CFLAGS)


hisq_invert_bench_f2_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_f2_CFLAGS = $(F2_CFLAGS)

hisq_invert_bench_d2_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_d2_CFLAGS = $(D2_CFLAGS)

hisq_invert_bench_f3_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_f3_CFLAGS = $(F3_CFLAGS)

hisq_invert_bench_d3_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_d3_CFLAGS = $(D3_CFLAGS)

hisq_invert_bench_fn_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_fn_CFLAGS = $(FN_CFLAGS)

hisq_invert_bench_dn_SOURCES = test_common.h test_common.c hisq_invert_bench.c
hisq_invert_bench_dn_CFLAGS = $(DN_CFLAGS)


hisq_force_bench_f2_SOURCES = test_common.h test_common.c hisq_force_bench.c
hisq_force_bench_f2_CFLAGS = $(F2_CFLAGS)

hisq_force_bench_d2_SOURCES = test_common.h test_common.c hisq_force_bench.c
hisq_force_bench_d2_CFLAGS = $(D2_CFLAGS)

hisq_force_bench_f3_SOURCES = test_common.h test_common.c hisq_force_bench.c
hisq_force_bench_f3_CFLAGS = $(F3_CFLAGS)

hisq_force_bench_d3_SOURCES = test_common.h test_common.c hisq_force_bench.c
hisq_force_bench_d3_CFLAGS = $(D3_CFLAGS)

hisq_force_bench_fn_SOURCES = test_common.h test_common.c hisq_force_bench.c
hisq_force_bench_fn_CFLAGS = $(FN_CFLAGS)

hisq_force_bench_dn_SOURCES = test_common.h test_common.c hisq_force_bench.c
hisq_force_bench_dn_CFLAGS = $(DN_CFLAGS)


wilson_test_f2_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_f2_CFLAGS = $(F2_CFLAGS)

wilson_test_d2_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_d2_CFLAGS = $(D2_CFLAGS)

wilson_test_f3_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_f3_CFLAGS = $(F3_CFLAGS)

wilson_test_d3_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_d3_CFLAGS = $(D3_CFLAGS)

wilson_test_fn_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_fn_CFLAGS = $(FN_CFLAGS)

wilson_test_dn_SOURCES = test_common.h test_common.c wilson_test.c
wilson_test_dn_CFLAGS = $(DN_CFLAGS)


wilson_invert_bench_f2_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_f2_CFLAGS = $(F2_CFLAGS)

wilson_invert_bench_d2_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_d2_CFLAGS = $(D2_CFLAGS)

wilson_invert_bench_f3_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_f3_CFLAGS = $(F3_CFLAGS)

wilson_invert_bench_d3_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_d3_CFLAGS = $(D3_CFLAGS)

wilson_invert_bench_fn_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_fn_CFLAGS = $(FN_CFLAGS)

wilson_invert_bench_dn_SOURCES = test_common.h test_common.c wilson_invert_bench.c
wilson_invert_bench_dn_CFLAGS = $(DN_CFLAGS)


wilson_force_bench_f2_SOURCES = test_common.h test_common.c wilson_force_bench.c
wilson_force_bench_f2_CFLAGS = $(F2_CFLAGS)

wilson_force_bench_d2_SOURCES = test_common.h test_common.c wilson_force_bench.c
wilson_force_bench_d2_CFLAGS = $(D2_CFLAGS)

wilson_force_bench_f3_SOURCES = test_common.h test_common.c wilson_force_bench.c
wilson_force_bench_f3_CFLAGS = $(F3_CFLAGS)

wilson_force_bench_d3_SOURCES = test_common.h test_common.c wilson_force_bench.c
wilson_force_bench_d3_CFLAGS = $(D3_CFLAGS)

wilson_force_bench_fn_SOURCES = test_common.h test_common.c wilson_force_bench.c
wilson_force_bench_fn_CFLAGS = $(FN_CFLAGS)

wilson_force_bench_dn_SOURCES = test_common.h test_common.c wilson_force_bench.c
wilson_force_bench_dn_CFLAGS = $(DN_CFLAGS)


dw_invert_bench_f2_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_f2_CFLAGS = $(F2_CFLAGS)

dw_invert_bench_d2_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_d2_CFLAGS = $(D2_CFLAGS)

dw_invert_bench_f3_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_f3_CFLAGS = $(F3_CFLAGS)

dw_invert_bench_d3_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_d3_CFLAGS = $(D3_CFLAGS)

dw_invert_bench_fn_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_fn_CFLAGS = $(FN_CFLAGS)

dw_invert_bench_dn_SOURCES = test_common.h test_common.c dw_invert_bench.c
dw_invert_bench_dn_CFLAGS = $(DN_CFLAGS)


test_reunit_utils_f3_SOURCES = test_common.h test_common.c test_reunit_utils.c
test_reunit_utils_f3_CFLAGS = $(F3_CFLAGS)

test_reunit_utils_d3_SOURCES = test_common.h test_common.c test_reunit_utils.c
test_reunit_utils_d3_CFLAGS = $(D3_CFLAGS)
