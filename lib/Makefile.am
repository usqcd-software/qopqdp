AM_CFLAGS = -I@top_srcdir@/include

lib_LIBRARIES = libqopqdp.a

libqopqdp_a_SOURCES = \
  common.c \
  layout_user.c \
  wilson_dslash.c \
  wilson_invert.c \
  mg/mg.c

NC1SOURCES = \
  pcstubs/common_f1.c pcstubs/common_d1.c \
  pcstubs/staples_f1.c pcstubs/staples_d1.c \
  pcstubs/gauge_deriv_f1.c pcstubs/gauge_deriv_d1.c \
  pcstubs/wilson_dslash_f1.c pcstubs/wilson_dslash_d1.c \
  pcstubs/wilson_invert_f1.c pcstubs/wilson_invert_d1.c \
  pcstubs/wilson_force_f1.c pcstubs/wilson_force_d1.c \
  pcstubs/wilson_mg_f1.c pcstubs/wilson_mg_d1.c \
  pcstubs/wilson_ops_f1.c pcstubs/wilson_ops_d1.c \
  pcstubs/invert_cg_V_f1.c pcstubs/invert_cg_V_d1.c \
  pcstubs/invert_cg_D_f1.c pcstubs/invert_cg_D_d1.c \
  pcstubs/invert_cg_vD_f1.c pcstubs/invert_cg_vD_d1.c \
  pcstubs/invert_bicgstab_D_f1.c pcstubs/invert_bicgstab_D_d1.c \
  pcstubs/invert_eigcg_V_f1.c pcstubs/invert_eigcg_V_d1.c \
  pcstubs/invert_eigcg_D_f1.c pcstubs/invert_eigcg_D_d1.c \
  pcstubs/invert_gmres2_D_f1.c pcstubs/invert_gmres2_D_d1.c \
  pcstubs/invert_utilities_V_f1.c pcstubs/invert_utilities_V_d1.c \
  pcstubs/invert_utilities_D_f1.c pcstubs/invert_utilities_D_d1.c \
  pcstubs/invert_utilities_vD_f1.c pcstubs/invert_utilities_vD_d1.c

NC2SOURCES = \
  pcstubs/common_f2.c pcstubs/common_d2.c \
  pcstubs/staples_f2.c pcstubs/staples_d2.c \
  pcstubs/gauge_deriv_f2.c pcstubs/gauge_deriv_d2.c \
  pcstubs/wilson_dslash_f2.c pcstubs/wilson_dslash_d2.c \
  pcstubs/wilson_invert_f2.c pcstubs/wilson_invert_d2.c \
  pcstubs/wilson_force_f2.c pcstubs/wilson_force_d2.c \
  pcstubs/wilson_mg_f2.c pcstubs/wilson_mg_d2.c \
  pcstubs/wilson_ops_f2.c pcstubs/wilson_ops_d2.c \
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
  pcstubs/staples_f3.c pcstubs/staples_d3.c \
  pcstubs/gauge_deriv_d3.c pcstubs/gauge_deriv_f3.c \
  pcstubs/wilson_dslash_f3.c pcstubs/wilson_dslash_d3.c \
  pcstubs/wilson_ifla_invert_f3.c pcstubs/wilson_ifla_invert_d3.c \
  pcstubs/wilson_invert_f3.c pcstubs/wilson_invert_d3.c \
  pcstubs/wilson_force_f3.c pcstubs/wilson_force_d3.c \
  pcstubs/wilson_mg_f3.c pcstubs/wilson_mg_d3.c \
  pcstubs/wilson_ops_f3.c pcstubs/wilson_ops_d3.c \
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

NCNSOURCES = \
  pcstubs/common_fn.c pcstubs/common_dn.c \
  pcstubs/staples_fn.c pcstubs/staples_dn.c \
  pcstubs/gauge_deriv_dn.c pcstubs/gauge_deriv_fn.c \
  pcstubs/wilson_dslash_fn.c pcstubs/wilson_dslash_dn.c \
  pcstubs/wilson_invert_fn.c pcstubs/wilson_invert_dn.c \
  pcstubs/wilson_force_fn.c pcstubs/wilson_force_dn.c \
  pcstubs/wilson_mg_fn.c pcstubs/wilson_mg_dn.c \
  pcstubs/wilson_ops_fn.c pcstubs/wilson_ops_dn.c \
  pcstubs/invert_cg_V_fn.c pcstubs/invert_cg_V_dn.c \
  pcstubs/invert_cg_D_fn.c pcstubs/invert_cg_D_dn.c \
  pcstubs/invert_cg_vD_fn.c pcstubs/invert_cg_vD_dn.c \
  pcstubs/invert_bicgstab_D_fn.c pcstubs/invert_bicgstab_D_dn.c \
  pcstubs/invert_eigcg_V_fn.c pcstubs/invert_eigcg_V_dn.c \
  pcstubs/invert_eigcg_D_fn.c pcstubs/invert_eigcg_D_dn.c \
  pcstubs/invert_gmres2_D_fn.c pcstubs/invert_gmres2_D_dn.c \
  pcstubs/invert_bicgstab_fn.c pcstubs/invert_bicgstab_dn.c \
  pcstubs/invert_cgls_fn.c pcstubs/invert_cgls_dn.c \
  pcstubs/invert_gcr_fn.c pcstubs/invert_gcr_dn.c \
  pcstubs/invert_utilities_V_fn.c pcstubs/invert_utilities_V_dn.c \
  pcstubs/invert_utilities_D_fn.c pcstubs/invert_utilities_D_dn.c \
  pcstubs/invert_utilities_vD_fn.c pcstubs/invert_utilities_vD_dn.c \
  pcstubs/mg_fn.c pcstubs/mg_dn.c \
  pcstubs/mg_dslash_fn.c pcstubs/mg_dslash_dn.c \
  pcstubs/vcycle_fn.c pcstubs/vcycle_dn.c

EXTRA_DIST = \
  common_p.c \
  staples_p.c \
  gauge_deriv_p.c \
  wilson_dslash_p.c \
  wilson_ifla_invert_p.c \
  wilson_invert_p.c \
  wilson_force_p.c \
  invert_cg_p.c \
  invert_cg_V_p.c \
  invert_cg_D_p.c \
  invert_cg_vD_p.c \
  invert_cg_blas_p.c \
  invert_bicgstab_p.c \
  invert_bicgstab_D_p.c \
  invert_eigcg_p.c \
  invert_eigcg_V_p.c \
  invert_eigcg_D_p.c \
  invert_gmres2_p.c \
  invert_gmres2_D_p.c \
  invert_utilities_p.c \
  invert_utilities_V_p.c \
  invert_utilities_D_p.c \
  invert_utilities_vD_p.c \
  linalg.h \
  mg

if NC1
  libqopqdp_a_SOURCES += $(NC1SOURCES)
endif
if NC2
  libqopqdp_a_SOURCES += $(NC2SOURCES)
endif
if NC3
  libqopqdp_a_SOURCES += $(NC3SOURCES)
endif
if NCN
  libqopqdp_a_SOURCES += $(NCNSOURCES)
endif
