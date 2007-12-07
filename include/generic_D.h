#define QOPPCV(x)		QOPPC(x ## _D)
#define Vector			QDP_DiracFermion
#define vIndexDef
#define create_V(r)		r = QDP_create_D()
#define destroy_V		QDP_destroy_D

#define V_eq_zero		QDP_D_eq_zero
#define V_eq_V			QDP_D_eq_D
#define V_peq_V			QDP_D_peq_D
#define V_eq_r_times_V		QDP_D_eq_r_times_D
#define V_peq_r_times_V		QDP_D_peq_r_times_D
#define V_meq_r_times_V		QDP_D_meq_r_times_D
#define V_peq_c_times_V		QDP_D_peq_c_times_D
#define V_meq_c_times_V		QDP_D_meq_c_times_D
#define V_eq_V_plus_V		QDP_D_eq_D_plus_D
#define V_eq_V_minus_V		QDP_D_eq_D_minus_D
#define V_eq_r_times_V_plus_V	QDP_D_eq_r_times_D_plus_D
#define V_eq_c_times_V_plus_V	QDP_D_eq_c_times_D_plus_D
#define r_eq_norm2_V		QDP_r_eq_norm2_D
#define r_eq_re_V_dot_V		QDP_r_eq_re_D_dot_D
#define c_eq_V_dot_V		QDP_c_eq_D_dot_D

#define V_veq_zero		QDP_D_veq_zero
#define V_veq_V			QDP_D_veq_D
#define V_vpeq_r_times_V	QDP_D_vpeq_r_times_D
#define V_vmeq_r_times_V	QDP_D_vmeq_r_times_D
#define V_vpeq_c_times_V	QDP_D_vpeq_c_times_D
#define V_veq_V_plus_V		QDP_D_veq_D_plus_D
#define V_veq_V_minus_V		QDP_D_veq_D_minus_D
#define V_veq_r_times_V_plus_V	QDP_D_veq_r_times_D_plus_D
#define V_veq_c_times_V_plus_V	QDP_D_veq_c_times_D_plus_D
#define r_veq_norm2_V		QDP_r_veq_norm2_D
#define r_veq_re_V_dot_V	QDP_r_veq_re_D_dot_D
