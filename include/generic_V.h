#define QOPPCV(x)		QOPPC(x ## _V)
#define Vector			QDP_ColorVector
#define vIndexDef
#define create_V(r)		r = QDP_create_V()
#define destroy_V		QDP_destroy_V

#define V_eq_zero		QDP_V_eq_zero
#define V_eq_V			QDP_V_eq_V
#define V_eq_r_times_V		QDP_V_eq_r_times_V
#define V_peq_r_times_V		QDP_V_peq_r_times_V
#define V_meq_r_times_V		QDP_V_meq_r_times_V
#define V_eq_V_plus_V		QDP_V_eq_V_plus_V
#define V_eq_V_minus_V		QDP_V_eq_V_minus_V
#define V_eq_r_times_V_plus_V	QDP_V_eq_r_times_V_plus_V
#define r_eq_norm2_V		QDP_r_eq_norm2_V
#define r_eq_re_V_dot_V		QDP_r_eq_re_V_dot_V

#define V_veq_zero		QDP_V_veq_zero
#define V_veq_V			QDP_V_veq_V
#define V_vpeq_r_times_V	QDP_V_vpeq_r_times_V
#define V_vmeq_r_times_V	QDP_V_vmeq_r_times_V
#define V_veq_V_plus_V		QDP_V_veq_V_plus_V
#define V_veq_V_minus_V		QDP_V_veq_V_minus_V
#define V_veq_r_times_V_plus_V	QDP_V_veq_r_times_V_plus_V
#define r_veq_norm2_V		QDP_r_veq_norm2_V
#define r_veq_re_V_dot_V	QDP_r_veq_re_V_dot_V
