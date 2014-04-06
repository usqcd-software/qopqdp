// DO NOT EDIT
// generated from qop_p_internal.h
#ifndef _QOP_D_INTERNAL_H
#define _QOP_D_INTERNAL_H

#ifdef HAVE_NCN
//#include <qop_mg_internal.h>

// BiCGStab

typedef struct {
  QDP_DN_ColorVector **r;
  QDP_DN_ColorVector **p;
  QDP_DN_ColorVector **v;
  QDP_DN_ColorVector **t;
  QDP_DN_ColorVector **Mx;
  int nv;
  int nc;
  int verbose;
  int indent;
} QOP_D_Bicgstab ;

typedef struct {
  void (*op)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  int itmax;
  QOP_D_Bicgstab *bicgstab;
  QDP_Subset sub;
  QDP_DN_ColorVector **outeo;
  QDP_DN_ColorVector **ineo;
  void (*project)(QDP_DN_ColorVector **ineo, QDP_DN_ColorVector **in, void *args);
  void (*reconstruct)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **outeo,
                      QDP_DN_ColorVector **in, void *args);
} QOP_D_BicgstabSolveArgs ;

QOP_D_Bicgstab * QOP_D_bicgstabInit (QDP_Lattice *lat, int nv, int nc);
void QOP_D_bicgstabFree (QOP_D_Bicgstab *bcg);
void QOP_D_bicgstabSet (QOP_D_Bicgstab *bcg, char *s, double v);
int QOP_D_bicgstabSolve (QOP_D_Bicgstab *bcg, QDP_DN_ColorVector *x[], QDP_DN_ColorVector *b[],
			 void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
			 void *Aargs,
			 void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
			 void *Margs,
			 double res, int itnlim, QDP_Subset sub);
int QOP_D_bicgstabSolveS (QOP_D_Bicgstab *bcg, QDP_DN_ColorVector *x[], QDP_DN_ColorVector *b[],
			  void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
			  void *Aargs,
			  void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
			  void *Margs,
			  int sign, double res, int itnlim, QDP_Subset sub);
int QOP_D_bicgstabSolveA (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
int QOP_D_bicgstabSolveEo (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);

// CGLS

typedef struct {
  QDP_DN_ColorVector **r1;
  QDP_DN_ColorVector **r2;
  QDP_DN_ColorVector **z;
  QDP_DN_ColorVector **p;
  QDP_DN_ColorVector **Ap;
  QDP_DN_ColorVector **Adr2;
  int nv;
  int nc;
  int verbose;
  int indent;
} QOP_D_Cgls ;

typedef struct {
  void (*op)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  double delta;
  int itmax;
  QOP_D_Cgls *cgls;
  QDP_Subset sub1;
  QDP_Subset sub2;
  QDP_DN_ColorVector **outeo;
  QDP_DN_ColorVector **ineo;
  void (*project)(QDP_DN_ColorVector **ineo, QDP_DN_ColorVector **in, void *args);
  void (*reconstruct)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **outeo,
                      QDP_DN_ColorVector **in, void *args);
} QOP_D_CglsSolveArgs ;

QOP_D_Cgls * QOP_D_cglsInit (QDP_Lattice *lat, int nv, int nc);
void QOP_D_cglsFree (QOP_D_Cgls *cgls);
void QOP_D_cglsSet (QOP_D_Cgls *cgls, char *s, double v);
int QOP_D_cglsSolve (QOP_D_Cgls *cgls, QDP_DN_ColorVector *x[],
		     QDP_DN_ColorVector *b1[], QDP_DN_ColorVector *b2[],
		     void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		     void *Aargs,
		     void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		     void *Margs,
		     double delta, double res1, double res2, int itnlim,
		     QDP_Subset sub1, QDP_Subset sub2);
int QOP_D_cgSolve (QOP_D_Cgls *cgls, QDP_DN_ColorVector *x[], QDP_DN_ColorVector *b[],
		   void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		   void *Aargs,
		   void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		   void *Margs,
		   double res, int itnlim, QDP_Subset sub);
int QOP_D_cglsSolveA1 (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
int QOP_D_cglsSolveA (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
int QOP_D_cglsSolveEo (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
int QOP_D_cglsSolveEo1 (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
int QOP_D_cgSolveA (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);

// GCR

typedef struct {
  QDP_DN_ColorVector **r;
  QDP_DN_ColorVector ***Mr;
  QDP_DN_ColorVector ***AMr;
  QDP_Lattice *lat;
  int nv;
  int nc;
  int ngcr;
  int restart;
  int verbose;
  int indent;
  int reuse;
  int nsmooth;
} QOP_D_Gcr;

typedef struct {
  void (*op)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*pop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *popargs;
  double res;
  int itmax;
  int ngcr;
  QOP_D_Gcr *gcr;
  QDP_Subset sub;
  QDP_DN_ColorVector **outeo;
  QDP_DN_ColorVector **ineo;
  void (*project)(QDP_DN_ColorVector **ineo, QDP_DN_ColorVector **in, void *args);
  void (*reconstruct)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **outeo,
                      QDP_DN_ColorVector **in, void *args);
} QOP_D_GcrSolveArgs;

QOP_D_Gcr *QOP_D_gcrInit (QDP_Lattice *lat, int nv, int nc, int ngcr);
void QOP_D_gcrFree (QOP_D_Gcr *gcr);
void QOP_D_gcrSet (QOP_D_Gcr *gcr, char *s, double v);
int QOP_D_gcrSolve (QOP_D_Gcr *gcr, QDP_DN_ColorVector *x[], QDP_DN_ColorVector *b[],
		    void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		    void *Aargs,
		    void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		    void *Margs,
		    double res, int itnlim, int ncgr, QDP_Subset subset);
int QOP_D_gcrSolve2 (QOP_D_Gcr *gcr, QDP_DN_ColorVector *x[], QDP_DN_ColorVector *b[],
		     void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], QOP_evenodd_t eo, void *args),
		     void *Aargs,
		     void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		     void *Margs,
		     double res, int itnlim, int ngcr);
int QOP_D_gcrlsSolve (QOP_D_Gcr *gcr, QDP_DN_ColorVector *x[],
		      QDP_DN_ColorVector *b1[], QDP_DN_ColorVector *b2[],
		      void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		      void *Aargs,
		      void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		      void *Margs,
		      double delta, double res1, double res2, int itnlim, int ngcr, QDP_Subset sub1, QDP_Subset sub2);
int QOP_D_gcrcgSolve (QOP_D_Gcr *gcr, QDP_DN_ColorVector *x[], QDP_DN_ColorVector *b[],
		      void Aop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		      void *Aargs,
		      void Mop(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
		      void *Margs,
		      double res, int itnlim, int ncgr, QDP_Subset subset);
int QOP_D_gcrSolveA (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
int QOP_D_gcrSolveEo (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
int QOP_D_gcrcgSolveA (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);

// V-cycle

typedef struct QOP_D_MgVcycleArgs {
  double s;
  double tpre;
  double tcoarse;
  double tpost;
  double delta;
  void (*op)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*cop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *copargs;
  void (*sop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *sopargs;
  int nv;
  int npre;
  int npost;
  int indent;
  int verbose;
  int count;
  QDP_Subset sub;
  QDP_Subset sub2;
  QOP_D_Gcr *gcr;
  QOP_D_Cgls *cgls;
  QDP_DN_ColorVector **r;
  QDP_DN_ColorVector **p;
  QDP_DN_ColorVector **Ap;
} QOP_D_MgVcycleArgs;

void QOP_D_mgVcycle(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);

#endif // HAVE_NCN

#if QOP_Precision == 'D'
#  include <qop_d_internal_generic.h>
#endif

#endif /* _QOP_D_INTERNAL_H */
