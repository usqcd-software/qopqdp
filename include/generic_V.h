//#define QOPPCV(x)		QOPPC(x ## _V)
#define QOPPCV(x)		QOP_ ## x ## _V
#define Vector			QDP_ColorVector
/* We need a Vector that is not pointerized in generic_vD.h */
#define Vector2			QDP_ColorVector
#define vector			QLA_ColorVector
#define _N                      1
#define vIndexDef
#define get_lattice_V(v)        QDP_get_lattice_V(v)
#define create_V(r)		r = QDP_create_V_L(lat)
#define destroy_V		QDP_destroy_V
#define expose_V                QDP_expose_V
#define reset_V                 QDP_reset_V

#define V_eq_funci              QDP_V_eq_funci
#define V_eq_funcit             QDP_V_eq_funcit

#define relnorm2_V(r,o,s)       QOP_relnorm2_V(&r, &o, s, 1);

#define insert_packed_V(r,a,s)  QDP_insert_packed_V(r,(QLA_ColorVector*)(a),s)
#define extract_packed_V(r,a,s) QDP_extract_packed_V((QLA_ColorVector*)(r),a,s)
#define csize_V                 QLA_Nc
#define first_qdp_object(x)     (x)

#define V_eq_zero		QDP_V_eq_zero
#define V_eq_V			QDP_V_eq_V
#define V_peq_V			QDP_V_peq_V
#define V_meq_V			QDP_V_meq_V
#define V_eq_V_plus_V		QDP_V_eq_V_plus_V
#define V_eq_V_minus_V		QDP_V_eq_V_minus_V

#define V_veq_zero		QDP_V_veq_zero
#define V_veq_V			QDP_V_veq_V
#define V_veq_V_plus_V		QDP_V_veq_V_plus_V
#define V_veq_V_minus_V		QDP_V_veq_V_minus_V

// scalars can be either single or double precision

#define r_eq_norm2_v(r,a)    { QLA_Real _r; QLA_R_eq_norm2_V(&_r,a); *(r) = _r; }
#define r_eq_norm2_V(r,a,s)  { QLA_Real _r; QDP_r_eq_norm2_V(&_r,a,s); *(r) = _r; }
#define r_eq_re_V_dot_V(r,a,b,s)  { QLA_Real _r; QDP_r_eq_re_V_dot_V(&_r,a,b,s); *(r) = _r; }
#define c_eq_V_dot_V(r,a,b,s)  { QLA_Complex _r; QDP_c_eq_V_dot_V(&_r,a,b,s); QLA_c_eq_c(*(r),_r); }

#define V_eq_r_times_V(r,a,b,s)	  { QLA_Real _a = *(a); QDP_V_eq_r_times_V(r,&_a,b,s); }
#define V_peq_r_times_V(r,a,b,s)  { QLA_Real _a = *(a); QDP_V_peq_r_times_V(r,&_a,b,s); }
#define V_meq_r_times_V(r,a,b,s)  { QLA_Real _a = *(a); QDP_V_meq_r_times_V(r,&_a,b,s); }
#define V_eq_c_times_V(r,a,b,s)  { QLA_Complex _a; QLA_c_eq_c(_a,*(a)); QDP_V_eq_c_times_V(r,&_a,b,s); }
#define V_peq_c_times_V(r,a,b,s)  { QLA_Complex _a; QLA_c_eq_c(_a,*(a)); QDP_V_peq_c_times_V(r,&_a,b,s); }
#define V_meq_c_times_V(r,a,b,s)  { QLA_Complex _a; QLA_c_eq_c(_a,*(a)); QDP_V_meq_c_times_V(r,&_a,b,s); }
#define V_eq_r_times_V_plus_V(r,a,b,c,s)  { QLA_Real _a = *(a); QDP_V_eq_r_times_V_plus_V(r,&_a,b,c,s); }
#define V_eq_c_times_V_plus_V(r,a,b,c,s)  { QLA_Complex _a; QLA_c_eq_c(_a,*(a)); QDP_V_eq_c_times_V_plus_V(r,&_a,b,c,s); }

// for vectorized routines, scalar arrays must be of the correct precision

#define r_veq_norm2_V     QDP_r_veq_norm2_V
#define r_veq_re_V_dot_V  QDP_r_veq_re_V_dot_V
#define c_veq_V_dot_V     QDP_c_veq_V_dot_V

#define V_vpeq_r_times_V  QDP_V_vpeq_r_times_V
#define V_vmeq_r_times_V  QDP_V_vmeq_r_times_V
#define V_vpeq_c_times_V  QDP_V_vpeq_c_times_V
#define V_veq_r_times_V_plus_V  QDP_V_veq_r_times_V_plus_V
#define V_veq_c_times_V_plus_V  QDP_V_veq_c_times_V_plus_V
