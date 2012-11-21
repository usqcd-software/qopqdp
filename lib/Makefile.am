AM_CFLAGS = -I@top_srcdir@/include

lib_LIBRARIES = libqopqdp.a

libqopqdp_a_SOURCES = \
  common.c \
  layout_user.c \
  asqtad_dslash.c \
  asqtad_invert.c \
  asqtad_force.c \
  hisq_links.c \
  hisq_force.c \
  wilson_dslash.c \
  wilson_invert.c \
  dw_dslash.c \
  dw_invert.c

NC2SOURCES = \
  pcstubs/common_f2.c pcstubs/common_d2.c \
  pcstubs/gauge_deriv_f2.c pcstubs/gauge_deriv_d2.c \
  pcstubs/symanzik_1loop_gauge_action_f2.c pcstubs/symanzik_1loop_gauge_action_d2.c \
  pcstubs/symanzik_1loop_gauge_force_f2.c pcstubs/symanzik_1loop_gauge_force_d2.c \
  pcstubs/asqtad_dslash_f2.c pcstubs/asqtad_dslash_d2.c \
  pcstubs/asqtad_invert_f2.c pcstubs/asqtad_invert_d2.c \
  pcstubs/asqtad_force_f2.c pcstubs/asqtad_force_d2.c \
  pcstubs/asqtad_force_asvec_f2.c pcstubs/asqtad_force_asvec_d2.c \
  pcstubs/asqtad_force_fnmat_f2.c pcstubs/asqtad_force_fnmat_d2.c \
  pcstubs/hisq_links_f2.c pcstubs/hisq_links_d2.c \
  pcstubs/hisq_reunit_f2.c pcstubs/hisq_reunit_d2.c \
  pcstubs/hisq_force_f2.c pcstubs/hisq_force_d2.c \
  pcstubs/hisq_force_fnmat_f2.c pcstubs/hisq_force_fnmat_d2.c \
  pcstubs/wilson_dslash_f2.c pcstubs/wilson_dslash_d2.c \
  pcstubs/wilson_invert_f2.c pcstubs/wilson_invert_d2.c \
  pcstubs/wilson_force_f2.c pcstubs/wilson_force_d2.c \
  pcstubs/dw_dslash_f2.c pcstubs/dw_dslash_d2.c \
  pcstubs/dw_invert_f2.c pcstubs/dw_invert_d2.c \
  pcstubs/invert_cg_V_f2.c pcstubs/invert_cg_V_d2.c \
  pcstubs/invert_cg_D_f2.c pcstubs/invert_cg_D_d2.c \
  pcstubs/invert_cg_vD_f2.c pcstubs/invert_cg_vD_d2.c \
  pcstubs/invert_bicgstab_D_f2.c pcstubs/invert_bicgstab_D_d2.c \
  pcstubs/invert_eigcg_V_f2.c pcstubs/invert_eigcg_V_d2.c \
  pcstubs/invert_eigcg_D_f2.c pcstubs/invert_eigcg_D_d2.c \
  pcstubs/invert_gmres2_D_f2.c pcstubs/invert_gmres2_D_d2.c \
  pcstubs/invert_utilities_V_f2.c pcstubs/invert_utilities_V_d2.c \
  pcstubs/invert_utilities_D_f2.c pcstubs/invert_utilities_D_d2.c \
  pcstubs/invert_utilities_vD_f2.c pcstubs/invert_utilities_vD_d2.c

NC3SOURCES = \
  pcstubs/common_f3.c pcstubs/common_d3.c \
  pcstubs/gauge_deriv_d3.c pcstubs/gauge_deriv_f3.c \
  pcstubs/symanzik_1loop_gauge_action_f3.c pcstubs/symanzik_1loop_gauge_action_d3.c \
  pcstubs/symanzik_1loop_gauge_force_f3.c pcstubs/symanzik_1loop_gauge_force_d3.c \
  pcstubs/asqtad_dslash_f3.c pcstubs/asqtad_dslash_d3.c \
  pcstubs/asqtad_invert_f3.c pcstubs/asqtad_invert_d3.c \
  pcstubs/asqtad_force_f3.c pcstubs/asqtad_force_d3.c \
  pcstubs/asqtad_force_asvec_f3.c pcstubs/asqtad_force_asvec_d3.c \
  pcstubs/asqtad_force_fnmat_f3.c pcstubs/asqtad_force_fnmat_d3.c \
  pcstubs/hisq_links_f3.c pcstubs/hisq_links_d3.c \
  pcstubs/hisq_reunit_f3.c pcstubs/hisq_reunit_d3.c \
  pcstubs/hisq_reunit_utilities_f3.c pcstubs/hisq_reunit_utilities_d3.c \
  pcstubs/hisq_force_f3.c pcstubs/hisq_force_d3.c \
  pcstubs/hisq_force_fnmat_f3.c pcstubs/hisq_force_fnmat_d3.c \
  pcstubs/wilson_dslash_f3.c pcstubs/wilson_dslash_d3.c \
  pcstubs/wilson_ifla_invert_f3.c pcstubs/wilson_ifla_invert_d3.c \
  pcstubs/wilson_invert_f3.c pcstubs/wilson_invert_d3.c \
  pcstubs/wilson_force_f3.c pcstubs/wilson_force_d3.c \
  pcstubs/dw_dslash_f3.c pcstubs/dw_dslash_d3.c \
  pcstubs/dw_invert_f3.c pcstubs/dw_invert_d3.c \
  pcstubs/invert_cg_V_f3.c pcstubs/invert_cg_V_d3.c \
  pcstubs/invert_cg_D_f3.c pcstubs/invert_cg_D_d3.c \
  pcstubs/invert_cg_vD_f3.c pcstubs/invert_cg_vD_d3.c \
  pcstubs/invert_bicgstab_D_f3.c pcstubs/invert_bicgstab_D_d3.c \
  pcstubs/invert_eigcg_V_f3.c pcstubs/invert_eigcg_V_d3.c \
  pcstubs/invert_eigcg_D_f3.c pcstubs/invert_eigcg_D_d3.c \
  pcstubs/invert_gmres2_D_f3.c pcstubs/invert_gmres2_D_d3.c \
  pcstubs/invert_utilities_V_f3.c pcstubs/invert_utilities_V_d3.c \
  pcstubs/invert_utilities_D_f3.c pcstubs/invert_utilities_D_d3.c \
  pcstubs/invert_utilities_vD_f3.c pcstubs/invert_utilities_vD_d3.c

EXTRA_DIST = \
  common_p.c \
  gauge_deriv_p.c \
  symanzik_1loop_gauge_action_p.c \
  symanzik_1loop_gauge_force2_p.c \
  asqtad_dslash_p.c \
  asqtad_invert_p.c \
  asqtad_force_p.c \
  asqtad_force_asvec_p.c \
  asqtad_force_fnmat3_p.c \
  hisq_force_fnmat2_p.c \
  hisq_force_p.c \
  hisq_links_p.c \
  hisq_reunit_p.c \
  hisq_reunit_p2.c \
  hisq_reunit_utilities_p.c \
  wilson_dslash_p.c \
  wilson_ifla_invert_p.c \
  wilson_invert_p.c \
  wilson_force_p.c \
  dw_dslash_p.c \
  dw_invert_p.c \
  invert_cg_p.c \
  invert_cg_blas_p.c \
  invert_bicgstab_p.c \
  invert_eigcg_p.c \
  invert_gmres2_p.c \
  invert_utilities_p.c \
  linalg.h

if NC2
  libqopqdp_a_SOURCES += $(NC2SOURCES)
endif
if NC3
  libqopqdp_a_SOURCES += $(NC3SOURCES)
endif
