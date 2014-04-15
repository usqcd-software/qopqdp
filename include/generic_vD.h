//#define QOPPCV(x)		QOPPC(x ## _vD)
#define QOPPCV(x)		QOP_ ## x ## _vD
//#define Vector			QDP_DiracFermion *
typedef QDP_DiracFermion * Vector;
#define Vector2			QDP_DiracFermion
#define _N                      _n
#define vIndexDef               , int _n
#define create_V(r)		{ r = malloc(_n*sizeof(Vector)); for(int _i=0; _i<_n; _i++) r[_i] = QDP_create_D(); }
#define destroy_V(r)		{ for(int _i=0; _i<_n; _i++) QDP_destroy_D(r[_i]); free(r); }
#define insert_packed_V(r,a,s)  { int _s=QDP_subset_len(s)*4*QLA_Nc; for(int _i=0; _i<_n; _i++) QDP_insert_packed_D((r)[_i],(void*)((a)+_s*_i),s); }
#define extract_packed_V(r,a,s)  { int _s=QDP_subset_len(s)*4*QLA_Nc; for(int _i=0; _i<_n; _i++) QDP_extract_packed_D((void*)((r)+_s*_i),(a)[_i],s); }
#define csize_V                 (4*QLA_Nc*_n)
#define first_qdp_object(x) ((x)[0])

/* Macros that are not vectorized */
#define vector                 QLA_DiracFermion
#define expose_V               QDP_expose_D
#define reset_V                QDP_reset_D
#define r_eq_norm2_v	       QLA_R_eq_norm2_D
#define V_eq_funci             QDP_D_eq_funci

#define relnorm2_V(r,o,s)      QOP_relnorm2_D(r, o, s, _n)

#define V_eq_zero(r,s)		{ for(int _i=0; _i<_n; _i++) QDP_D_eq_zero(r[_i],s); }
#define V_eq_V(r,a,s)		QDP_D_veq_D(r,a,s,_n)
#define V_meq_V(r,a,s)		QDP_D_vmeq_D(r,a,s,_n)
#define V_eq_V_plus_V(r,a,b,s)	QDP_D_veq_D_plus_D(r,a,b,s,_n)
#define V_eq_V_minus_V(r,a,b,s)	QDP_D_veq_D_minus_D(r,a,b,s,_n)

// all scalars can be either single or double precision

#define r_eq_norm2_V(r,a,s)  { QLA_Real _r[_n]; QDP_r_veq_norm2_D(_r,a,s,_n); *(r) = 0; for(int _i=0; _i<_n; _i++) *(r) += _r[_i]; }
#define r_eq_re_V_dot_V(r,a,b,s)  { QLA_Real _r[_n]; QDP_r_veq_re_D_dot_D(_r,a,b,s,_n); *(r) = 0; for(int _i=0; _i<_n; _i++) *(r) += _r[_i]; }
#define c_eq_V_dot_V(r,a,b,s)  { QLA_Complex _r[_n]; QDP_c_veq_D_dot_D(_r,a,b,s,_n); QLA_c_eq_r(*(r),0); for(int _i=0; _i<_n; _i++) QLA_c_peq_c(*(r),_r[_i]);}

#define V_eq_r_times_V(r,a,b,s)  { QLA_Real _a[_n]; for(int _i=0; _i<_n; _i++) _a[_i] = *(a); QDP_D_veq_r_times_D(r,_a,b,s,_n); }
#define V_peq_r_times_V(r,a,b,s)  { QLA_Real _a[_n]; for(int _i=0; _i<_n; _i++) _a[_i] = *(a); QDP_D_vpeq_r_times_D(r,_a,b,s,_n); }
#define V_meq_r_times_V(r,a,b,s)  { QLA_Real _a[_n]; for(int _i=0; _i<_n; _i++) _a[_i] = *(a); QDP_D_vmeq_r_times_D(r,_a,b,s,_n); }
#define V_eq_c_times_V(r,a,b,s)  { QLA_Complex _a[_n]; for(int _i=0; _i<_n; _i++) QLA_c_eq_c(_a[_i],*(a)); QDP_D_veq_c_times_D(r,_a,b,s,_n); }
#define V_peq_c_times_V(r,a,b,s)  { QLA_Complex _a[_n]; for(int _i=0; _i<_n; _i++) QLA_c_eq_c(_a[_i],*(a)); QDP_D_vpeq_c_times_D(r,_a,b,s,_n); }
#define V_eq_r_times_V_plus_V(r,a,b,c,s)  { QLA_Real _a[_n]; for(int _i=0; _i<_n; _i++) _a[_i] = *(a); QDP_D_veq_r_times_D_plus_D(r,_a,b,c,s,_n); }
#define V_eq_c_times_V_plus_V(r,a,b,c,s)  { QLA_Complex _a[_n]; for(int _i=0; _i<_n; _i++) QLA_c_eq_c(_a[_i],*(a)); QDP_D_veq_c_times_D_plus_D(r,_a,b,c,s,_n); }
