#define QOPPCV(x)		QOPPC(x ## _D)
#define Vector			QDP_DiracFermion
/* We need a Vector that is not pointerized in generic_vD.h */
#define Vector2			QDP_DiracFermion
#define vector			QLA_DiracFermion
#define _N                      1
#define vIndexDef
#define create_V(r)		r = QDP_create_D()
#define destroy_V		QDP_destroy_D
#define expose_V                QDP_expose_D
#define reset_V                 QDP_reset_D

#define V_eq_funci              QDP_D_eq_funci

#define relnorm2_V(r,o,s)       QOPPCV(relnorm2)(&r, &o, s, 1);

#define insert_packed_V(r,a,s)  QDP_insert_packed_D(r,(QLA_DiracFermion*)(a),s)
#define extract_packed_V(r,a,s) QDP_extract_packed_D((QLA_DiracFermion*)(r),a,s)
#define csize_V                 12

#define V_eq_zero		QDP_D_eq_zero
#define V_eq_V			QDP_D_eq_D
#define V_peq_V			QDP_D_peq_D
#define V_meq_V			QDP_D_meq_D
#define V_eq_V_plus_V		QDP_D_eq_D_plus_D
#define V_eq_V_minus_V		QDP_D_eq_D_minus_D

#define V_veq_zero		QDP_D_veq_zero
#define V_veq_V			QDP_D_veq_D
#define V_veq_V_plus_V		QDP_D_veq_D_plus_D
#define V_veq_V_minus_V		QDP_D_veq_D_minus_D

// scalars can be either single or double precision

#define r_eq_norm2_v(r,a)    { QLA_Real _r; QLA_R_eq_norm2_D(&_r,a); *(r) = _r; }
#define r_eq_norm2_V(r,a,s)  { QLA_Real _r; QDP_r_eq_norm2_D(&_r,a,s); *(r) = _r; }
#define r_eq_re_V_dot_V(r,a,b,s)  { QLA_Real _r; QDP_r_eq_re_D_dot_D(&_r,a,b,s); *(r) = _r; }
#define c_eq_V_dot_V(r,a,b,s)  { QLA_Complex _r; QDP_c_eq_D_dot_D(&_r,a,b,s); QLA_c_eq_c(*(r),_r); }

#define V_eq_r_times_V(r,a,b,s)	  { QLA_Real _a = *(a); QDP_D_eq_r_times_D(r,&_a,b,s); }
#define V_eq_c_times_V(r,a,b,s)   { QLA_Complex _a; QLA_c_eq_c(_a,*(a)); QDP_D_eq_c_times_D(r,&_a,b,s); }
#define V_peq_r_times_V(r,a,b,s)  { QLA_Real _a = *(a); QDP_D_peq_r_times_D(r,&_a,b,s); }
#define V_meq_r_times_V(r,a,b,s)  { QLA_Real _a = *(a); QDP_D_meq_r_times_D(r,&_a,b,s); }
#define V_peq_c_times_V(r,a,b,s)  { QLA_Complex _a; QLA_c_eq_c(_a,*(a)); QDP_D_peq_c_times_D(r,&_a,b,s); }
#define V_meq_c_times_V(r,a,b,s)  { QLA_Complex _a; QLA_c_eq_c(_a,*(a)); QDP_D_meq_c_times_D(r,&_a,b,s); }
#define V_eq_r_times_V_plus_V(r,a,b,c,s)  { QLA_Real _a = *(a); QDP_D_eq_r_times_D_plus_D(r,&_a,b,c,s); }

// for vectorized routines, scalar arrays must be of the correct precision

#define r_veq_norm2_V     QDP_r_veq_norm2_D
#define r_veq_re_V_dot_V  QDP_r_veq_re_D_dot_D
#define c_veq_V_dot_V     QDP_c_veq_D_dot_D

#define V_vpeq_r_times_V  QDP_D_vpeq_r_times_D
#define V_vmeq_r_times_V  QDP_D_vmeq_r_times_D
#define V_vpeq_c_times_V  QDP_D_vpeq_c_times_D
#define V_veq_r_times_V_plus_V  QDP_D_veq_r_times_D_plus_D
#define V_veq_c_times_V_plus_V  QDP_D_veq_c_times_D_plus_D
