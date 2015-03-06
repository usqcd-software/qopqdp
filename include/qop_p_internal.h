#ifndef _QOP_P_INTERNAL_H
#define _QOP_P_INTERNAL_H

#ifdef HAVE_NCN
//#include <qop_mg_internal.h>

// BiCGStab

typedef struct {
  QDP_N_ColorVector **r;
  QDP_N_ColorVector **p;
  QDP_N_ColorVector **v;
  QDP_N_ColorVector **t;
  QDP_N_ColorVector **Mx;
  int nv;
  int nc;
  int verbose;
  int indent;
} QOP_P_Bicgstab ;

typedef struct {
  void (*op)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  int itmax;
  QOP_P_Bicgstab *bicgstab;
  QDP_Subset sub;
  QDP_N_ColorVector **outeo;
  QDP_N_ColorVector **ineo;
  void (*project)(QDP_N_ColorVector **ineo, QDP_N_ColorVector **in, void *args);
  void (*reconstruct)(QDP_N_ColorVector **out, QDP_N_ColorVector **outeo,
                      QDP_N_ColorVector **in, void *args);
} QOP_P_BicgstabSolveArgs ;

QOP_P_Bicgstab * QOP_P_bicgstabInit (QDP_Lattice *lat, int nv, int nc);
void QOP_P_bicgstabFree (QOP_P_Bicgstab *bcg);
void QOP_P_bicgstabSet (QOP_P_Bicgstab *bcg, char *s, double v);
int QOP_P_bicgstabSolve (QOP_P_Bicgstab *bcg, QDP_N_ColorVector *x[], QDP_N_ColorVector *b[],
			 void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
			 void *Aargs,
			 void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
			 void *Margs,
			 double res, int itnlim, QDP_Subset sub);
int QOP_P_bicgstabSolveS (QOP_P_Bicgstab *bcg, QDP_N_ColorVector *x[], QDP_N_ColorVector *b[],
			  void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
			  void *Aargs,
			  void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
			  void *Margs,
			  int sign, double res, int itnlim, QDP_Subset sub);
int QOP_P_bicgstabSolveA (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
int QOP_P_bicgstabSolveEo (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);

// CGLS

typedef struct {
  QDP_N_ColorVector **r1;
  QDP_N_ColorVector **r2;
  QDP_N_ColorVector **z;
  QDP_N_ColorVector **p;
  QDP_N_ColorVector **Ap;
  QDP_N_ColorVector **Adr2;
  int nv;
  int nc;
  int verbose;
  int indent;
} QOP_P_Cgls ;

typedef struct {
  void (*op)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  double delta;
  int itmax;
  QOP_P_Cgls *cgls;
  QDP_Subset sub1;
  QDP_Subset sub2;
  QDP_N_ColorVector **outeo;
  QDP_N_ColorVector **ineo;
  void (*project)(QDP_N_ColorVector **ineo, QDP_N_ColorVector **in, void *args);
  void (*reconstruct)(QDP_N_ColorVector **out, QDP_N_ColorVector **outeo,
                      QDP_N_ColorVector **in, void *args);
} QOP_P_CglsSolveArgs ;

QOP_P_Cgls * QOP_P_cglsInit (QDP_Lattice *lat, int nv, int nc);
void QOP_P_cglsFree (QOP_P_Cgls *cgls);
void QOP_P_cglsSet (QOP_P_Cgls *cgls, char *s, double v);
int QOP_P_cglsSolve (QOP_P_Cgls *cgls, QDP_N_ColorVector *x[],
		     QDP_N_ColorVector *b1[], QDP_N_ColorVector *b2[],
		     void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		     void *Aargs,
		     void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		     void *Margs,
		     double delta, double res1, double res2, int itnlim,
		     QDP_Subset sub1, QDP_Subset sub2);
int QOP_P_cgSolve (QOP_P_Cgls *cgls, QDP_N_ColorVector *x[], QDP_N_ColorVector *b[],
		   void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		   void *Aargs,
		   void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		   void *Margs,
		   double res, int itnlim, QDP_Subset sub);
int QOP_P_cglsSolveA1 (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
int QOP_P_cglsSolveA (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
int QOP_P_cglsSolveEo (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
int QOP_P_cglsSolveEo1 (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
int QOP_P_cgSolveA (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);

// GCR

typedef struct {
  QDP_N_ColorVector **r;
  QDP_N_ColorVector ***Mr;
  QDP_N_ColorVector ***AMr;
  QDP_Lattice *lat;
  int nv;
  int nc;
  int ngcr;
  int restart;
  int verbose;
  int indent;
  int reuse;
  int nsmooth;
} QOP_IP_Gcr;

typedef struct {
  void (*op)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  int itmax;
  int ngcr;
  QOP_P_Gcr *gcr;
  QDP_Subset sub;
  QDP_N_ColorVector **outeo;
  QDP_N_ColorVector **ineo;
  void (*project)(QDP_N_ColorVector **ineo, QDP_N_ColorVector **in, void *args);
  void (*reconstruct)(QDP_N_ColorVector **out, QDP_N_ColorVector **outeo,
                      QDP_N_ColorVector **in, void *args);
} QOP_IP_GcrSolveArgs;

QOP_P_Gcr *QOP_P_gcrInit (QDP_Lattice *lat, int nv, int nc, int ngcr);
void QOP_P_gcrFree (QOP_P_Gcr *gcr);
void QOP_P_gcrSet (QOP_P_Gcr *gcr, char *s, double v);
int QOP_P_gcrSolve (QOP_P_Gcr *gcr, QDP_N_ColorVector *x[], QDP_N_ColorVector *b[],
		    void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		    void *Aargs,
		    void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		    void *Margs,
		    double res, int itnlim, int ncgr, QDP_Subset subset);
int QOP_P_gcrSolve2 (QOP_P_Gcr *gcr, QDP_N_ColorVector *x[], QDP_N_ColorVector *b[],
		     void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], QOP_evenodd_t eo, void *args),
		     void *Aargs,
		     void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		     void *Margs,
		     double res, int itnlim, int ngcr);
int QOP_P_gcrlsSolve (QOP_P_Gcr *gcr, QDP_N_ColorVector *x[],
		      QDP_N_ColorVector *b1[], QDP_N_ColorVector *b2[],
		      void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		      void *Aargs,
		      void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		      void *Margs,
		      double delta, double res1, double res2, int itnlim, int ngcr, QDP_Subset sub1, QDP_Subset sub2);
int QOP_P_gcrcgSolve (QOP_P_Gcr *gcr, QDP_N_ColorVector *x[], QDP_N_ColorVector *b[],
		      void Aop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		      void *Aargs,
		      void Mop(QDP_N_ColorVector *Ax[], QDP_N_ColorVector *x[], int sign, void *args),
		      void *Margs,
		      double res, int itnlim, int ncgr, QDP_Subset subset);
int QOP_P_gcrSolveA (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
int QOP_P_gcrSolveEo (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
int QOP_P_gcrcgSolveA (QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);

// V-cycle

typedef struct QOP_IP_MgVcycleArgs {
  double s;
  double tpre;
  double tcoarse;
  double tpost;
  double delta;
  void (*op)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*cop)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *copargs;
  void (*sop)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *sopargs;
  int nv;
  int npre;
  int npost;
  int indent;
  int verbose;
  int count;
  QDP_Subset sub;
  QDP_Subset sub2;
  QOP_P_Gcr *gcr;
  QOP_P_Cgls *cgls;
  QDP_N_ColorVector **r;
  QDP_N_ColorVector **p;
  QDP_N_ColorVector **Ap;
} QOP_IP_MgVcycleArgs;

void QOP_IP_mgVcycle(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);

#endif // HAVE_NCN

#if QOP_Precision == _QOP_Precision
#  include <qop_p_internal_generic.h>
#endif

#ifdef HAVE_QLL

void IP(setup_qll)(QDP_Lattice *lat);
void * IP(get_qll_layout)(void);
void IP(toQDP)(QLA_Real *xx, QDP_Lattice *lat, QLA_Real *yy, void *l, int nelem);
void IP(fromQDP)(QLA_Real *yy, void *l, QLA_Real *xx, QDP_Lattice *lat, int nelem);

#endif // HAVE_QLL

#endif /* _QOP_P_INTERNAL_H */
