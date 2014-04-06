// DO NOT EDIT
// generated from qop_p_internal.h
#ifndef _QOP_F_INTERNAL_H
#define _QOP_F_INTERNAL_H

#ifdef HAVE_NCN
//#include <qop_mg_internal.h>

// BiCGStab

typedef struct {
  QDP_FN_ColorVector **r;
  QDP_FN_ColorVector **p;
  QDP_FN_ColorVector **v;
  QDP_FN_ColorVector **t;
  QDP_FN_ColorVector **Mx;
  int nv;
  int nc;
  int verbose;
  int indent;
} QOP_F_Bicgstab ;

typedef struct {
  void (*op)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  int itmax;
  QOP_F_Bicgstab *bicgstab;
  QDP_Subset sub;
  QDP_FN_ColorVector **outeo;
  QDP_FN_ColorVector **ineo;
  void (*project)(QDP_FN_ColorVector **ineo, QDP_FN_ColorVector **in, void *args);
  void (*reconstruct)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **outeo,
                      QDP_FN_ColorVector **in, void *args);
} QOP_F_BicgstabSolveArgs ;

QOP_F_Bicgstab * QOP_F_bicgstabInit (QDP_Lattice *lat, int nv, int nc);
void QOP_F_bicgstabFree (QOP_F_Bicgstab *bcg);
void QOP_F_bicgstabSet (QOP_F_Bicgstab *bcg, char *s, double v);
int QOP_F_bicgstabSolve (QOP_F_Bicgstab *bcg, QDP_FN_ColorVector *x[], QDP_FN_ColorVector *b[],
			 void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
			 void *Aargs,
			 void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
			 void *Margs,
			 double res, int itnlim, QDP_Subset sub);
int QOP_F_bicgstabSolveS (QOP_F_Bicgstab *bcg, QDP_FN_ColorVector *x[], QDP_FN_ColorVector *b[],
			  void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
			  void *Aargs,
			  void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
			  void *Margs,
			  int sign, double res, int itnlim, QDP_Subset sub);
int QOP_F_bicgstabSolveA (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
int QOP_F_bicgstabSolveEo (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

// CGLS

typedef struct {
  QDP_FN_ColorVector **r1;
  QDP_FN_ColorVector **r2;
  QDP_FN_ColorVector **z;
  QDP_FN_ColorVector **p;
  QDP_FN_ColorVector **Ap;
  QDP_FN_ColorVector **Adr2;
  int nv;
  int nc;
  int verbose;
  int indent;
} QOP_F_Cgls ;

typedef struct {
  void (*op)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  double delta;
  int itmax;
  QOP_F_Cgls *cgls;
  QDP_Subset sub1;
  QDP_Subset sub2;
  QDP_FN_ColorVector **outeo;
  QDP_FN_ColorVector **ineo;
  void (*project)(QDP_FN_ColorVector **ineo, QDP_FN_ColorVector **in, void *args);
  void (*reconstruct)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **outeo,
                      QDP_FN_ColorVector **in, void *args);
} QOP_F_CglsSolveArgs ;

QOP_F_Cgls * QOP_F_cglsInit (QDP_Lattice *lat, int nv, int nc);
void QOP_F_cglsFree (QOP_F_Cgls *cgls);
void QOP_F_cglsSet (QOP_F_Cgls *cgls, char *s, double v);
int QOP_F_cglsSolve (QOP_F_Cgls *cgls, QDP_FN_ColorVector *x[],
		     QDP_FN_ColorVector *b1[], QDP_FN_ColorVector *b2[],
		     void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		     void *Aargs,
		     void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		     void *Margs,
		     double delta, double res1, double res2, int itnlim,
		     QDP_Subset sub1, QDP_Subset sub2);
int QOP_F_cgSolve (QOP_F_Cgls *cgls, QDP_FN_ColorVector *x[], QDP_FN_ColorVector *b[],
		   void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		   void *Aargs,
		   void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		   void *Margs,
		   double res, int itnlim, QDP_Subset sub);
int QOP_F_cglsSolveA1 (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
int QOP_F_cglsSolveA (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
int QOP_F_cglsSolveEo (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
int QOP_F_cglsSolveEo1 (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
int QOP_F_cgSolveA (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

// GCR

typedef struct {
  QDP_FN_ColorVector **r;
  QDP_FN_ColorVector ***Mr;
  QDP_FN_ColorVector ***AMr;
  QDP_Lattice *lat;
  int nv;
  int nc;
  int ngcr;
  int restart;
  int verbose;
  int indent;
  int reuse;
  int nsmooth;
} QOP_F_Gcr;

typedef struct {
  void (*op)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  int itmax;
  int ngcr;
  QOP_F_Gcr *gcr;
  QDP_Subset sub;
  QDP_FN_ColorVector **outeo;
  QDP_FN_ColorVector **ineo;
  void (*project)(QDP_FN_ColorVector **ineo, QDP_FN_ColorVector **in, void *args);
  void (*reconstruct)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **outeo,
                      QDP_FN_ColorVector **in, void *args);
} QOP_F_GcrSolveArgs;

QOP_F_Gcr *QOP_F_gcrInit (QDP_Lattice *lat, int nv, int nc, int ngcr);
void QOP_F_gcrFree (QOP_F_Gcr *gcr);
void QOP_F_gcrSet (QOP_F_Gcr *gcr, char *s, double v);
int QOP_F_gcrSolve (QOP_F_Gcr *gcr, QDP_FN_ColorVector *x[], QDP_FN_ColorVector *b[],
		    void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		    void *Aargs,
		    void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		    void *Margs,
		    double res, int itnlim, int ncgr, QDP_Subset subset);
int QOP_F_gcrSolve2 (QOP_F_Gcr *gcr, QDP_FN_ColorVector *x[], QDP_FN_ColorVector *b[],
		     void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], QOP_evenodd_t eo, void *args),
		     void *Aargs,
		     void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		     void *Margs,
		     double res, int itnlim, int ngcr);
int QOP_F_gcrlsSolve (QOP_F_Gcr *gcr, QDP_FN_ColorVector *x[],
		      QDP_FN_ColorVector *b1[], QDP_FN_ColorVector *b2[],
		      void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		      void *Aargs,
		      void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		      void *Margs,
		      double delta, double res1, double res2, int itnlim, int ngcr, QDP_Subset sub1, QDP_Subset sub2);
int QOP_F_gcrcgSolve (QOP_F_Gcr *gcr, QDP_FN_ColorVector *x[], QDP_FN_ColorVector *b[],
		      void Aop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		      void *Aargs,
		      void Mop(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		      void *Margs,
		      double res, int itnlim, int ncgr, QDP_Subset subset);
int QOP_F_gcrSolveA (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
int QOP_F_gcrSolveEo (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
int QOP_F_gcrcgSolveA (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

// V-cycle

typedef struct QOP_F_MgVcycleArgs {
  double s;
  double tpre;
  double tcoarse;
  double tpost;
  double delta;
  void (*op)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*cop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *copargs;
  void (*sop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *sopargs;
  int nv;
  int npre;
  int npost;
  int indent;
  int verbose;
  int count;
  QDP_Subset sub;
  QDP_Subset sub2;
  QOP_F_Gcr *gcr;
  QOP_F_Cgls *cgls;
  QDP_FN_ColorVector **r;
  QDP_FN_ColorVector **p;
  QDP_FN_ColorVector **Ap;
} QOP_F_MgVcycleArgs;

void QOP_F_mgVcycle(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

#endif // HAVE_NCN

#if QOP_Precision == 'F'
#  include <qop_f_internal_generic.h>
#endif

#endif /* _QOP_F_INTERNAL_H */
