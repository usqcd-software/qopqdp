AM_CFLAGS = -I@top_srcdir@/include

lib_LIBRARIES = libqopqdp.a

libqopqdp_a_SOURCES = \
  common.c common_f.c common_d.c \
  symanzik_1loop_gauge_force_f.c symanzik_1loop_gauge_force_d.c \
  asqtad_dslash_f.c asqtad_dslash_d.c \
  asqtad_invert.c asqtad_invert_f.c asqtad_invert_d.c \
  asqtad_force.c asqtad_force_f.c asqtad_force_d.c \
  asqtad_force_fnmat_f.c asqtad_force_fnmat_d.c \
  hisq_links.c hisq_links_f.c hisq_links_d.c \
  hisq_reunit_f.c hisq_reunit_d.c \
  hisq_reunit_utilities_f.c hisq_reunit_utilities_d.c\
  hisq_force_f.c hisq_force_d.c hisq_force.c \
  hisq_force_fnmat_f.c hisq_force_fnmat_d.c \
  wilson_dslash.c wilson_dslash_f.c wilson_dslash_d.c \
  wilson_invert.c wilson_invert_f.c wilson_invert_d.c \
  dw_dslash.c dw_dslash_f.c dw_dslash_d.c \
  dw_invert.c dw_invert_f.c dw_invert_d.c \
  invert_cg_V_f.c invert_cg_V_d.c \
  invert_cg_D_f.c invert_cg_D_d.c \
  invert_cg_vD_f.c invert_cg_vD_d.c \
  invert_bicgstab_D_f.c invert_bicgstab_D_d.c \
  invert_eigcg_V_f.c invert_eigcg_V_d.c \
  invert_eigcg_D_f.c invert_eigcg_D_d.c \
  invert_gmres2_D_f.c invert_gmres2_D_d.c \
  invert_utilities_V_d.c invert_utilities_V_f.c \
  invert_utilities_D_d.c invert_utilities_D_f.c \
  invert_utilities_vD_d.c invert_utilities_vD_f.c

common_f.o: common_p.c
common_d.o: common_p.c
symanzik_1loop_gauge_force_f.o: symanzik_1loop_gauge_force_p.c
symanzik_1loop_gauge_force_d.o: symanzik_1loop_gauge_force_p.c
asqtad_dslash_f.o: asqtad_dslash_p.c
asqtad_dslash_d.o: asqtad_dslash_p.c
asqtad_invert_f.o: asqtad_invert_p.c
asqtad_invert_d.o: asqtad_invert_p.c
asqtad_force_f.o: asqtad_force_p.c
asqtad_force_d.o: asqtad_force_p.c
asqtad_force_fnmat_f.o: asqtad_force_fnmat_p.c
asqtad_force_fnmat_d.o: asqtad_force_fnmat_p.c
hisq_links_f.o: hisq_links_p.c
hisq_links_d.o: hisq_links_p.c
hisq_reunit_f.o: hisq_reunit_p.c
hisq_reunit_d.o: hisq_reunit_p.c
hisq_reunit_utilities_f.o: hisq_reunit_utilities_p.c
hisq_reunit_utilities_d.o: hisq_reunit_utilities_p.c
hisq_force_f.o: hisq_force_p.c
hisq_force_d.o: hisq_force_p.c
hisq_force_fnmat_f.o: hisq_force_fnmat_p.c
hisq_force_fnmat_d.o: hisq_force_fnmat_p.c
wilson_dslash_f.o: wilson_dslash_p.c
wilson_dslash_d.o: wilson_dslash_p.c
wilson_invert_f.o: wilson_invert_p.c
wilson_invert_d.o: wilson_invert_p.c
dw_dslash_f.o: dw_dslash_p.c
dw_dslash_d.o: dw_dslash_p.c
dw_invert_f.o: dw_invert_p.c
dw_invert_d.o: dw_invert_p.c
invert_cg_V_f.o: invert_cg_p.c
invert_cg_V_d.o: invert_cg_p.c
invert_cg_D_f.o: invert_cg_p.c
invert_cg_D_d.o: invert_cg_p.c
invert_cg_vD_f.o: invert_cg_p.c
invert_cg_vD_d.o: invert_cg_p.c
invert_bicgstab_D_f.o: invert_bicgstab_p.c
invert_bicgstab_D_d.o: invert_bicgstab_p.c
invert_eigcg_V_f.o: invert_eigcg_p.c
invert_eigcg_V_d.o: invert_eigcg_p.c
invert_eigcg_D_f.o: invert_eigcg_p.c
invert_eigcg_D_d.o: invert_eigcg_p.c
invert_gmres2_D_f.o: invert_gmres2_p.c
invert_gmres2_D_d.o: invert_gmres2_p.c
invert_utilities_V_f.o: invert_utilities_p.c
invert_utilities_V_d.o: invert_utilities_p.c
invert_utilities_D_f.o: invert_utilities_p.c
invert_utilities_D_d.o: invert_utilities_p.c
invert_utilities_vD_f.o: invert_utilities_p.c
invert_utilities_vD_d.o: invert_utilities_p.c

EXTRA_DIST = \
  common_p.c \
  symanzik_1loop_gauge_force_p.c \
  asqtad_dslash_p.c \
  asqtad_invert_p.c \
  asqtad_force_p.c \
  asqtad_force_fnmat_p.c \
  hisq_force_fnmat_p.c \
  hisq_force_p.c \
  hisq_links_p.c \
  hisq_reunit_p.c \
  hisq_reunit_utilities_p.c \
  wilson_dslash_p.c \
  wilson_invert_p.c \
  dw_dslash_p.c \
  dw_invert_p.c \
  invert_cg_p.c \
  invert_cg_blas_p.c \
  invert_bicgstab_p.c \
  invert_eigcg_p.c \
  invert_gmres2_p.c \
  invert_utilities_p.c \
  linalg.h
