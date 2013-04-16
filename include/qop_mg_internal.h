typedef struct {
  QOP_F3_FermionLinksWilson *wil;
  QLA_F_Real kappa;
} QOP_F3_WilArgs;

typedef struct {
  QOP_D3_FermionLinksWilson *wil;
  QLA_D_Real kappa;
} QOP_D3_WilArgs;

void QOP_F3_wilsonDslash (QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *in,
			  QOP_F3_FermionLinksWilson *wil, QLA_F_Real kappa, int sign,
			  QOP_evenodd_t pout, QOP_evenodd_t pin);
void QOP_F3_wilsonDiaginv (QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *in,
			   QOP_F3_FermionLinksWilson *wil, QLA_F_Real kappa,
			   QOP_evenodd_t pout);
void QOP_F3_wilsonDslashEO (QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *in,
			    QOP_F3_FermionLinksWilson *wil, QLA_F_Real kappa, int sign, QOP_evenodd_t par);
void QOP_F3_wilsonDslashEOS (QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *in,
			     QOP_F3_FermionLinksWilson *wil, QLA_F_Real kappa, int sign, QOP_evenodd_t par);
void QOP_F3_wilsonDslashEOH (QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *in,
			     QOP_F3_FermionLinksWilson *wil, QLA_F_Real kappa, int sign, QOP_evenodd_t par);

void QOP_F3_V1eqD (QDP_FN_ColorVector *v[1], QDP_F3_DiracFermion *d, QDP_Subset sub);
void QOP_F3_DeqV1 (QDP_DiracFermion *d, QDP_FN_ColorVector *v[1], QDP_Subset sub);
void QOP_F3_V2eqD (QDP_FN_ColorVector *v[2], QDP_F3_DiracFermion *d, QDP_Subset sub);
void QOP_F3_DeqV2 (QDP_DiracFermion *d, QDP_FN_ColorVector *v[2], QDP_Subset sub);
void QOP_F3_wilDV1 (QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_F3_wilDV2 (QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F3_wilPV1 (QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_F3_wilPV2 (QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F3_wilPNEV2 (QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F3_wilEoV1 (QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_F3_wilEoV2 (QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F3_wilEoProjectD (QDP_F3_DiracFermion *ineo, QDP_F3_DiracFermion *in, QOP_F3_WilArgs *w);
void QOP_F3_wilEoReconstructD (QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *outeo,
			       QDP_F3_DiracFermion *in, QOP_F3_WilArgs *w);
void QOP_F3_wilEoProjectV1 (QDP_FN_ColorVector *ineo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_F3_wilEoReconstructV1 (QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *outeo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_F3_wilEoReconstructPV1 (QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *outeo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_F3_wilEoProjectV2 (QDP_FN_ColorVector *ineo[2], QDP_FN_ColorVector *in[2], void *args);
void QOP_F3_wilEoReconstructV2 (QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *outeo[2], QDP_FN_ColorVector *in[2], void *args);
void QOP_F3_wilEoReconstructPV2 (QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *outeo[2], QDP_FN_ColorVector *in[2], void *args);

void QOP_D3_wilsonDslash (QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in,
			  QOP_D3_FermionLinksWilson *wil, QLA_D_Real kappa, int sign,
			  QOP_evenodd_t pout, QOP_evenodd_t pin);
void QOP_D3_wilsonDiaginv (QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in,
			   QOP_D3_FermionLinksWilson *wil, QLA_D_Real kappa,
			   QOP_evenodd_t pout);
void QOP_D3_wilsonDslashEO (QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in,
			    QOP_D3_FermionLinksWilson *wil, QLA_D_Real kappa, int sign, QOP_evenodd_t par);
void QOP_D3_wilsonDslashEOS (QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in,
			     QOP_D3_FermionLinksWilson *wil, QLA_D_Real kappa, int sign, QOP_evenodd_t par);
void QOP_D3_wilsonDslashEOH (QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in,
			     QOP_D3_FermionLinksWilson *wil, QLA_D_Real kappa, int sign, QOP_evenodd_t par);

void QOP_D3_V1eqD (QDP_DN_ColorVector *v[1], QDP_D3_DiracFermion *d, QDP_Subset sub);
void QOP_D3_DeqV1 (QDP_DiracFermion *d, QDP_DN_ColorVector *v[1], QDP_Subset sub);
void QOP_D3_V2eqD (QDP_DN_ColorVector *v[2], QDP_D3_DiracFermion *d, QDP_Subset sub);
void QOP_D3_DeqV2 (QDP_DiracFermion *d, QDP_DN_ColorVector *v[2], QDP_Subset sub);
void QOP_D3_wilDV1 (QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_D3_wilDV2 (QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D3_wilPV1 (QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_D3_wilPV2 (QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D3_wilPNEV2 (QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D3_wilEoV1 (QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_D3_wilEoV2 (QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D3_wilEoProjectD (QDP_D3_DiracFermion *ineo, QDP_D3_DiracFermion *in, QOP_D3_WilArgs *w);
void QOP_D3_wilEoReconstructD (QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *outeo,
			       QDP_D3_DiracFermion *in, QOP_D3_WilArgs *w);
void QOP_D3_wilEoProjectV1 (QDP_DN_ColorVector *ineo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_D3_wilEoReconstructV1 (QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *outeo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_D3_wilEoReconstructPV1 (QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *outeo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_D3_wilEoProjectV2 (QDP_DN_ColorVector *ineo[2], QDP_DN_ColorVector *in[2], void *args);
void QOP_D3_wilEoReconstructV2 (QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *outeo[2], QDP_DN_ColorVector *in[2], void *args);
void QOP_D3_wilEoReconstructPV2 (QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *outeo[2], QDP_DN_ColorVector *in[2], void *args);

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
} QOP_D_Gcr ;

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
} QOP_D_GcrSolveArgs ;

QOP_D_Gcr * QOP_D_gcrInit (QDP_Lattice *lat, int nv, int nc, int ngcr);
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

// MG block

typedef struct {
  int *sites[2];
  int nsites[2];
} QOP_MgLocalBlock ;

typedef struct {
  double trestrict;
  double tprolong;
  QDP_Lattice *fine;
  QDP_Lattice *coarse;
  QDP_Subset *fs;
  QDP_Subset *cs;
  int ns;
  int local;
  int nlb;
  int rcount;
  int pcount;
  QOP_MgLocalBlock *lb;
} QOP_MgBlock ;

QOP_MgBlock * QOP_mgCreateBlock (QDP_Lattice *fine, int *block);
QOP_MgBlock * QOP_mgCreateBlockFromLattice (QDP_Lattice *fine, QDP_Lattice *coarse);
void QOP_mgFreeBlock (QOP_MgBlock *mgb);

// MG args

typedef struct {
  QOP_MgBlock *mgb;
  QDP_FN_ColorVector ***rv;
  QDP_FN_ColorVector ***pv;
  int cnc;
  int fnc;
  int nv;
  int maxnv;
  int we_malloced_rv;
  int we_malloced_pv;
} QOP_F_MgArgs ;

typedef void QOP_F_MgOp(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

typedef struct {
  QOP_F_MgArgs *mga;
  int (*op)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  //QOP_F_MgOp *op;
  void *opargs;
  QDP_FN_ColorVector **cin;
  QDP_FN_ColorVector **cout;
  QOP_evenodd_t fpar;
} QOP_F_MgF2cOpArgs ;

typedef struct {
  QOP_F_MgArgs *mga;
  QOP_F_MgOp *op;
  void *opargs;
  QDP_FN_ColorVector **fin;
  QDP_FN_ColorVector **fout;
  QOP_evenodd_t fpar;
} QOP_F_MgC2fOpArgs ;

QOP_F_MgArgs *QOP_F_mgCreateArgs(QOP_MgBlock *mgb, int cnc, int fnc, int nv, QDP_FN_ColorVector ***rv, QDP_FN_ColorVector ***pv);
void QOP_F_mgFreeArgs (QOP_F_MgArgs *mga);
void QOP_F_mgF2c (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, QOP_F_MgArgs *w, QOP_evenodd_t par);
void QOP_F_mgC2f (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, QOP_F_MgArgs *w, QOP_evenodd_t par);
void QOP_F_mgF2cOp (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
void QOP_F_mgC2fOp (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

void QOP_F_mgOrthoVn (QDP_FN_ColorVector ***vv, int nv, int i0, int n, QDP_Subset sub);
void QOP_F_mgOrtho (QDP_FN_ColorVector *cv[], int nv, QOP_MgBlock *mgb);
void QOP_F_mgOrthoVec (QDP_FN_ColorVector *cv[], int i, QOP_MgBlock *mgb, int norm);
void QOP_F_mgOrtho2 (QDP_FN_ColorVector *cv1[], QDP_FN_ColorVector *cv2[], int nv, QOP_MgBlock *mgb);
void QOP_F_rRitzHarm (QDP_FN_ColorVector *v[], int n, int nv, QOP_F_MgOp *op, void *opargs, QLA_F_Complex evs[], int sort, int norm, QDP_Subset sub);
void QOP_F_mgOrthoSort (QDP_FN_ColorVector *cv[], int imin, int n, QOP_MgBlock *mgb, double min[], double ave[], double max[]);

void QOP_F_mgRestrict (QDP_FN_ColorVector *cv[], QDP_FN_ColorVector *fv[], int nv,
		       QDP_FN_ColorVector *pv[], int cnc, int fnc, QOP_MgBlock *mgb, QOP_evenodd_t par);
void QOP_F_mgProlong (QDP_FN_ColorVector *fv[], QDP_FN_ColorVector *cv[], int nv,
		      QDP_FN_ColorVector *pv[], int cnc, int fnc, QOP_MgBlock *mgb, QOP_evenodd_t par);
void QOP_F_mgTestCoarse (QOP_F_MgArgs *mga, QDP_FN_ColorVector **vf);

typedef struct {
  QOP_MgBlock *mgb;
  QDP_DN_ColorVector ***rv;
  QDP_DN_ColorVector ***pv;
  int cnc;
  int fnc;
  int nv;
  int maxnv;
  int we_malloced_rv;
  int we_malloced_pv;
} QOP_D_MgArgs ;

typedef void QOP_D_MgOp (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int 
			 sign, void *args);

typedef struct {
  QOP_D_MgArgs *mga;
  int (*op)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  //QOP_D_MgOp *op;
  void *opargs;
  QDP_DN_ColorVector **cin;
  QDP_DN_ColorVector **cout;
  QOP_evenodd_t fpar;
} QOP_D_MgF2cOpArgs ;

typedef struct {
  QOP_D_MgArgs *mga;
  QOP_D_MgOp *op;
  void *opargs;
  QDP_DN_ColorVector **fin;
  QDP_DN_ColorVector **fout;
  QOP_evenodd_t fpar;
} QOP_D_MgC2fOpArgs ;

QOP_D_MgArgs * QOP_D_mgCreateArgs (QOP_MgBlock *mgb, int cnc, int fnc, int nv, QDP_DN_ColorVector ***rv, QDP_DN_ColorVector ***pv);
void QOP_D_mgFreeArgs (QOP_D_MgArgs *mga);
void QOP_D_mgF2c (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, QOP_D_MgArgs *w, QOP_evenodd_t par);
void QOP_D_mgC2f (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, QOP_D_MgArgs *w, QOP_evenodd_t par);
void QOP_D_mgF2cOp (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
void QOP_D_mgC2fOp (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);

void QOP_D_mgOrthoVn (QDP_DN_ColorVector ***vv, int nv, int i0, int n, QDP_Subset sub);
void QOP_D_mgOrtho (QDP_DN_ColorVector *cv[], int nv, QOP_MgBlock *mgb);
void QOP_D_mgOrthoVec (QDP_DN_ColorVector *cv[], int i, QOP_MgBlock *mgb, int norm);
void QOP_D_mgOrtho2 (QDP_DN_ColorVector *cv1[], QDP_DN_ColorVector *cv2[], int nv, QOP_MgBlock *mgb);
void QOP_D_rRitzHarm (QDP_DN_ColorVector *v[], int n, int nv, QOP_D_MgOp *op, void *opargs, QLA_D_Complex evs[], int sort, int norm, QDP_Subset sub);
void QOP_D_mgOrthoSort (QDP_DN_ColorVector *cv[], int imin, int n, QOP_MgBlock *mgb, double min[], double ave[], double max[]);

void QOP_D_mgRestrict (QDP_DN_ColorVector *cv[], QDP_DN_ColorVector *fv[], int nv,
		       QDP_DN_ColorVector *pv[], int cnc, int fnc, QOP_MgBlock *mgb, QOP_evenodd_t par);
void QOP_D_mgProlong (QDP_DN_ColorVector *fv[], QDP_DN_ColorVector *cv[], int nv,
		      QDP_DN_ColorVector *pv[], int cnc, int fnc, QOP_MgBlock *mgb, QOP_evenodd_t par);
void QOP_D_mgTestCoarse (QOP_D_MgArgs *mga, QDP_DN_ColorVector **vf);

// MG dslash args

typedef struct {
  QDP_FN_ColorVector ***temp;
  QDP_FN_ColorMatrix ****links;
  int nvout;
  int nvin;
  int nc;
  int nshifts;
  int npaths;
  int *paths;
  QDP_Lattice *lat;
  QDP_Subset even;
  QDP_Subset odd;
  QDP_Subset all;
  QDP_Shift *shifts;
  QDP_ShiftDir *fb;
} QOP_F_MgDslashArgs ;

QOP_F_MgDslashArgs * QOP_F_mgCreateDslash (int nvout, int nvin, int nc, int npaths, int *paths, QDP_Lattice *lat);
void QOP_F_mgFreeDslash (QOP_F_MgDslashArgs *da);
void QOP_F_mgCloneOp (void op(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
		      void *args, QOP_F_MgDslashArgs *da);
void QOP_F_mgSetDiaginv (QOP_F_MgDslashArgs *da);
void QOP_F_mgTestDslash (void op(QDP_FN_ColorVector *Ax[], QDP_FN_ColorVector *x[], int sign, void *args),
			 void *args, QOP_F_MgDslashArgs *da);
void QOP_F_mgDslash (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args,
		     QOP_evenodd_t parout, QOP_evenodd_t parin);
void QOP_F_mgDslashAll (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args);
void QOP_F_mgDiag (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args,
		   QOP_evenodd_t par);
void QOP_F_mgDiaginv (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args,
		      QOP_evenodd_t par);
void QOP_F_mgDslashP (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args,
		      QOP_evenodd_t parout, QOP_evenodd_t parin);
void QOP_F_mgDslashPAll (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args);
void QOP_F_mgDslashEo (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args);
void QOP_F_mgDslashEoProject (QDP_FN_ColorVector *ineo[], QDP_FN_ColorVector *in[], void *args);
void QOP_F_mgDslashEoReconstruct (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *outeo[],
				  QDP_FN_ColorVector *in[], void *args);
void QOP_F_mgDslashEoReconstructP (QDP_FN_ColorVector *out[], QDP_FN_ColorVector *outeo[],
				   QDP_FN_ColorVector *in[], void *args);

typedef struct {
  QDP_DN_ColorVector ***temp;
  QDP_DN_ColorMatrix ****links;
  int nvout;
  int nvin;
  int nc;
  int nshifts;
  int npaths;
  int *paths;
  QDP_Lattice *lat;
  QDP_Subset even;
  QDP_Subset odd;
  QDP_Subset all;
  QDP_Shift *shifts;
  QDP_ShiftDir *fb;
} QOP_D_MgDslashArgs;

QOP_D_MgDslashArgs * QOP_D_mgCreateDslash (int nvout, int nvin, int nc, int npaths, int *paths, QDP_Lattice *lat);
void QOP_D_mgFreeDslash (QOP_D_MgDslashArgs *da);
void QOP_D_mgCloneOp (void op(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[],
			      int sign, void *args),
		      void *args, QOP_D_MgDslashArgs *da);
void QOP_D_mgSetDiaginv (QOP_D_MgDslashArgs *da);
void QOP_D_mgTestDslash (void op(QDP_DN_ColorVector *Ax[], QDP_DN_ColorVector *x[], int sign, void *args),
			 void *args, QOP_D_MgDslashArgs *da);
void QOP_D_mgDslash (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *in[], int sign, void *args,
		     QOP_evenodd_t parout, QOP_evenodd_t parin);
void QOP_D_mgDslashAll (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *in[], int sign, void *args);
void QOP_D_mgDiag (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *in[], int sign, void *args,
		   QOP_evenodd_t par);
void QOP_D_mgDiaginv (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *in[], int sign, void *args,
		      QOP_evenodd_t par);
void QOP_D_mgDslashP (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *in[], int sign, void *args,
		      QOP_evenodd_t parout, QOP_evenodd_t parin);
void QOP_D_mgDslashPAll (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *in[], int sign, void *args);
void QOP_D_mgDslashEo (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *in[], int sign, void *args);
void QOP_D_mgDslashEoProject (QDP_DN_ColorVector *ineo[], QDP_DN_ColorVector *in[], void *args);
void QOP_D_mgDslashEoReconstruct (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *outeo[],
				  QDP_DN_ColorVector *in[], void *args);
void QOP_D_mgDslashEoReconstructP (QDP_DN_ColorVector *out[], QDP_DN_ColorVector *outeo[],
				   QDP_DN_ColorVector *in[], void *args);

// V-cycle

typedef struct {
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
} QOP_F_MgVcycleArgs ;

void QOP_F_mgVcycle (QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

void QOP_F_eignhVN (int n, int nv, QDP_FN_ColorVector *v[n][nv], int nc,
		    QDP_FN_ColorVector *v0[nv], QOP_F_MgOp *op, void *opargs,
		    QDP_Subset sub, int maxits);
int QOP_F_svdVN (int n, int nv, QDP_FN_ColorVector *v[n][nv], int nc,
		 QDP_FN_ColorVector *v0[nv], QOP_F_MgOp *op, void *opargs,
		 QDP_Subset sub, QDP_Subset sub2, int maxits);

typedef struct {
  double s;
  double tpre;
  double tcoarse;
  double tpost;
  double delta;
  void (*op)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *
	     args);
  void *opargs;
  void (*cop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void 
	      *args);
  void *copargs;
  void (*sop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void 
	      *args);
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
} QOP_D_MgVcycleArgs ;

void QOP_D_mgVcycle (QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);

void QOP_D_eignhVN (int n, int nv, QDP_DN_ColorVector *v[n][nv], int nc,
		    QDP_DN_ColorVector *v0[nv], QOP_D_MgOp *op, void *opargs,
		    QDP_Subset sub, int maxits);
int QOP_D_svdVN (int n, int nv, QDP_DN_ColorVector *v[n][nv], int nc,
		 QDP_DN_ColorVector *v0[nv], QOP_D_MgOp *op, void *opargs,
		 QDP_Subset sub, QDP_Subset sub2, int maxits);

// Wilson MG

typedef struct {
  // general parameters
  int ndim;
  int *lattice_size;
  int fnc;
  int nv;
  int nvecs;
  double setup_res;
  double setup_change_fac;
  int setup_maxit;
  int setup_nvecs;
  int verbose;
  // Vcycle parameters
  int npre;
  int npost;
  QLA_F_Real scale;
  QOP_F_MgOp *vcop;
  void *vcopargs;
  QOP_F_Gcr *gcrf;
  QOP_F_MgF2cOpArgs *fcoa; //
  QOP_F_MgVcycleArgs *vca; //
  // coarse operator variables
  QDP_Lattice *lattice;
  QOP_MgBlock *mgblock;
  QOP_F_MgArgs *mgargs;
  QOP_F_MgDslashArgs *dargs;
  QOP_F_MgC2fOpArgs *cfoa; //
  QOP_F_MgOp *nvop;
  void *nvopargs;
  QOP_F_MgOp *cfop;
  void *cfopargs;
  // coarse solve variables
  double cres;
  int itmax;
  int ngcr;
  QOP_F_Gcr *gcrc;
  QOP_F_GcrSolveArgs *sa; //
  QOP_F_MgOp *csop;
  void *csopargs;
  int created;
} QOP_WilMgLevel ;

struct QOP_WilsonMgStruct {
  QOP_F3_FermionLinksWilson *wilF;
  QOP_D3_FermionLinksWilson *wilD;
  QLA_F_Real kappa;
  QLA_F_Real kappanv;
  QOP_F3_WilArgs vcwaF;
  QOP_F3_WilArgs nvwaF;
  int nlevels;
  QOP_WilMgLevel *mg;
  int verbose;
  int profile;
  int itmax;
  QOP_F_Gcr *gcrF;
  QOP_D_Gcr *gcrD;
  int ngcr;
};


#if QOP_Precision == 'F'

#define QOP_WilArgs QOP_F3_WilArgs
#define QOP_wilsonNew QOP_F3_wilsonNew
#define QOP_wilsonFree QOP_F3_wilsonFree
#define QOP_wilsonSet QOP_F3_wilsonSet
#define QOP_wilsonSetCudaDevice QOP_F3_wilsonSetCudaDevice
#define QOP_wilsonLoadLinks QOP_F3_wilsonLoadLinks
#define QOP_wilsonDslash QOP_F3_wilsonDslash
#define QOP_wilsonDiaginv QOP_F3_wilsonDiaginv
#define QOP_wilsonDslashEO QOP_F3_wilsonDslashEO
#define QOP_wilsonDslashEOS QOP_F3_wilsonDslashEOS
#define QOP_wilsonDslashEOH QOP_F3_wilsonDslashEOH
#define QOP_wilsonDslashVec QOP_F3_wilsonDslashVec
#define QOP_wilsonDslashEOVec QOP_F3_wilsonDslashEOVec
#define QOP_wilsonDslashEOSVec QOP_F3_wilsonDslashEOSVec
#define QOP_wilsonDslashEOHVec QOP_F3_wilsonDslashEOHVec
#define QOP_wilsonInvert QOP_F3_wilsonInvert
#define QOP_wilsonEigen QOP_F3_wilsonEigen
#define QOP_V1eqD QOP_F3_V1eqD
#define QOP_DeqV1 QOP_F3_DeqV1
#define QOP_V2eqD QOP_F3_V2eqD
#define QOP_DeqV2 QOP_F3_DeqV2
#define QOP_wilDV1 QOP_F3_wilDV1
#define QOP_wilDV2 QOP_F3_wilDV2
#define QOP_wilPV1 QOP_F3_wilPV1
#define QOP_wilPV2 QOP_F3_wilPV2
#define QOP_wilPNEV2 QOP_F3_wilPNEV2
#define QOP_wilEoV1 QOP_F3_wilEoV1
#define QOP_wilEoV2 QOP_F3_wilEoV2
#define QOP_wilEoProjectD QOP_F3_wilEoProjectD
#define QOP_wilEoReconstructD QOP_F3_wilEoReconstructD
#define QOP_wilEoProjectV1 QOP_F3_wilEoProjectV1
#define QOP_wilEoReconstructV1 QOP_F3_wilEoReconstructV1
#define QOP_wilEoReconstructPV1 QOP_F3_wilEoReconstructPV1
#define QOP_wilEoProjectV2 QOP_F3_wilEoProjectV2
#define QOP_wilEoReconstructV2 QOP_F3_wilEoReconstructV2
#define QOP_wilEoReconstructPV2 QOP_F3_wilEoReconstructPV2
#define QOP_getDiagEOS QOP_F3_getDiagEOS

#define QOP_MgArgs QOP_F_MgArgs
#define QOP_MgOp QOP_F_MgOp
#define QOP_MgF2cOpArgs QOP_F_MgF2cOpArgs
#define QOP_MgC2fOpArgs QOP_F_MgC2fOpArgs
#define QOP_mgCreateArgs QOP_F_mgCreateArgs
#define QOP_mgFreeArgs QOP_F_mgFreeArgs
#define QOP_mgF2c QOP_F_mgF2c
#define QOP_mgC2f QOP_F_mgC2f
#define QOP_mgF2cOp QOP_F_mgF2cOp
#define QOP_mgC2fOp QOP_F_mgC2fOp
#define QOP_mgOrthoVn QOP_F_mgOrthoVn
#define QOP_mgOrtho QOP_F_mgOrtho
#define QOP_mgOrthoVec QOP_F_mgOrthoVec
#define QOP_mgOrtho2 QOP_F_mgOrtho2
#define QOP_rRitzHarm QOP_F_rRitzHarm
#define QOP_mgOrthoSort QOP_F_mgOrthoSort
#define QOP_mgRestrict QOP_F_mgRestrict
#define QOP_mgProlong QOP_F_mgProlong
#define QOP_mgTestCoarse QOP_F_mgTestCoarse
#define QOP_MgDslashArgs QOP_F_MgDslashArgs
#define QOP_mgCreateDslash QOP_F_mgCreateDslash
#define QOP_mgFreeDslash QOP_F_mgFreeDslash
#define QOP_mgCloneOp QOP_F_mgCloneOp
#define QOP_mgSetDiaginv QOP_F_mgSetDiaginv
#define QOP_mgTestDslash QOP_F_mgTestDslash
#define QOP_mgDslash QOP_F_mgDslash
#define QOP_mgDslashAll QOP_F_mgDslashAll
#define QOP_mgDiag QOP_F_mgDiag
#define QOP_mgDiaginv QOP_F_mgDiaginv
#define QOP_mgDslashP QOP_F_mgDslashP
#define QOP_mgDslashPAll QOP_F_mgDslashPAll
#define QOP_mgDslashEo QOP_F_mgDslashEo
#define QOP_mgDslashEoProject QOP_F_mgDslashEoProject
#define QOP_mgDslashEoReconstruct QOP_F_mgDslashEoReconstruct
#define QOP_mgDslashEoReconstructP QOP_F_mgDslashEoReconstructP
#define QOP_Gcr QOP_F_Gcr
#define QOP_GcrSolveArgs QOP_F_GcrSolveArgs
#define QOP_gcrInit QOP_F_gcrInit
#define QOP_gcrFree QOP_F_gcrFree
#define QOP_gcrSet QOP_F_gcrSet
#define QOP_gcrSolve QOP_F_gcrSolve
#define QOP_gcrSolve2 QOP_F_gcrSolve2
#define QOP_gcrlsSolve QOP_F_gcrlsSolve
#define QOP_gcrcgSolve QOP_F_gcrcgSolve
#define QOP_gcrSolveA QOP_F_gcrSolveA
#define QOP_gcrSolveEo QOP_F_gcrSolveEo
#define QOP_gcrcgSolveA QOP_F_gcrcgSolveA
#define QOP_Cgls QOP_F_Cgls
#define QOP_CglsSolveArgs QOP_F_CglsSolveArgs
#define QOP_cglsInit QOP_F_cglsInit
#define QOP_cglsFree QOP_F_cglsFree
#define QOP_cglsSet QOP_F_cglsSet
#define QOP_cglsSolve QOP_F_cglsSolve
#define QOP_cgSolve QOP_F_cgSolve
#define QOP_cglsSolveA1 QOP_F_cglsSolveA1
#define QOP_cglsSolveA QOP_F_cglsSolveA
#define QOP_cglsSolveEo QOP_F_cglsSolveEo
#define QOP_cglsSolveEo1 QOP_F_cglsSolveEo1
#define QOP_cgSolveA QOP_F_cgSolveA
#define QOP_Bicgstab QOP_F_Bicgstab
#define QOP_BicgstabSolveArgs QOP_F_BicgstabSolveArgs
#define QOP_bicgstabInit QOP_F_bicgstabInit
#define QOP_bicgstabFree QOP_F_bicgstabFree
#define QOP_bicgstabSet QOP_F_bicgstabSet
#define QOP_bicgstabSolve QOP_F_bicgstabSolve
#define QOP_bicgstabSolveS QOP_F_bicgstabSolveS
#define QOP_bicgstabSolveA QOP_F_bicgstabSolveA
#define QOP_bicgstabSolveEo QOP_F_bicgstabSolveEo
#define QOP_MgVcycleArgs QOP_F_MgVcycleArgs
#define QOP_mgVcycle QOP_F_mgVcycle
#define QOP_eignhVN QOP_F_eignhVN
#define QOP_svdVN QOP_F_svdVN

#else

#define QOP_WilArgs QOP_D3_WilArgs
#define QOP_wilsonNew QOP_D3_wilsonNew
#define QOP_wilsonFree QOP_D3_wilsonFree
#define QOP_wilsonSet QOP_D3_wilsonSet
#define QOP_wilsonSetCudaDevice QOP_D3_wilsonSetCudaDevice
#define QOP_wilsonLoadLinks QOP_D3_wilsonLoadLinks
#define QOP_wilsonDslash QOP_D3_wilsonDslash
#define QOP_wilsonDiaginv QOP_D3_wilsonDiaginv
#define QOP_wilsonDslashEO QOP_D3_wilsonDslashEO
#define QOP_wilsonDslashEOS QOP_D3_wilsonDslashEOS
#define QOP_wilsonDslashEOH QOP_D3_wilsonDslashEOH
#define QOP_wilsonDslashVec QOP_D3_wilsonDslashVec
#define QOP_wilsonDslashEOVec QOP_D3_wilsonDslashEOVec
#define QOP_wilsonDslashEOSVec QOP_D3_wilsonDslashEOSVec
#define QOP_wilsonDslashEOHVec QOP_D3_wilsonDslashEOHVec
#define QOP_wilsonInvert QOP_D3_wilsonInvert
#define QOP_wilsonEigen QOP_D3_wilsonEigen
#define QOP_V1eqD QOP_D3_V1eqD
#define QOP_DeqV1 QOP_D3_DeqV1
#define QOP_V2eqD QOP_D3_V2eqD
#define QOP_DeqV2 QOP_D3_DeqV2
#define QOP_wilDV1 QOP_D3_wilDV1
#define QOP_wilDV2 QOP_D3_wilDV2
#define QOP_wilPV1 QOP_D3_wilPV1
#define QOP_wilPV2 QOP_D3_wilPV2
#define QOP_wilPNEV2 QOP_D3_wilPNEV2
#define QOP_wilEoV1 QOP_D3_wilEoV1
#define QOP_wilEoV2 QOP_D3_wilEoV2
#define QOP_wilEoProjectD QOP_D3_wilEoProjectD
#define QOP_wilEoReconstructD QOP_D3_wilEoReconstructD
#define QOP_wilEoProjectV1 QOP_D3_wilEoProjectV1
#define QOP_wilEoReconstructV1 QOP_D3_wilEoReconstructV1
#define QOP_wilEoReconstructPV1 QOP_D3_wilEoReconstructPV1
#define QOP_wilEoProjectV2 QOP_D3_wilEoProjectV2
#define QOP_wilEoReconstructV2 QOP_D3_wilEoReconstructV2
#define QOP_wilEoReconstructPV2 QOP_D3_wilEoReconstructPV2
#define QOP_getDiagEOS QOP_D3_getDiagEOS

#define QOP_MgArgs QOP_D_MgArgs
#define QOP_MgOp QOP_D_MgOp
#define QOP_MgF2cOpArgs QOP_D_MgF2cOpArgs
#define QOP_MgC2fOpArgs QOP_D_MgC2fOpArgs
#define QOP_mgCreateArgs QOP_D_mgCreateArgs
#define QOP_mgFreeArgs QOP_D_mgFreeArgs
#define QOP_mgF2c QOP_D_mgF2c
#define QOP_mgC2f QOP_D_mgC2f
#define QOP_mgF2cOp QOP_D_mgF2cOp
#define QOP_mgC2fOp QOP_D_mgC2fOp
#define QOP_mgOrthoVn QOP_D_mgOrthoVn
#define QOP_mgOrtho QOP_D_mgOrtho
#define QOP_mgOrthoVec QOP_D_mgOrthoVec
#define QOP_mgOrtho2 QOP_D_mgOrtho2
#define QOP_rRitzHarm QOP_D_rRitzHarm
#define QOP_mgOrthoSort QOP_D_mgOrthoSort
#define QOP_mgRestrict QOP_D_mgRestrict
#define QOP_mgProlong QOP_D_mgProlong
#define QOP_mgTestCoarse QOP_D_mgTestCoarse
#define QOP_MgDslashArgs QOP_D_MgDslashArgs
#define QOP_mgCreateDslash QOP_D_mgCreateDslash
#define QOP_mgFreeDslash QOP_D_mgFreeDslash
#define QOP_mgCloneOp QOP_D_mgCloneOp
#define QOP_mgSetDiaginv QOP_D_mgSetDiaginv
#define QOP_mgTestDslash QOP_D_mgTestDslash
#define QOP_mgDslash QOP_D_mgDslash
#define QOP_mgDslashAll QOP_D_mgDslashAll
#define QOP_mgDiag QOP_D_mgDiag
#define QOP_mgDiaginv QOP_D_mgDiaginv
#define QOP_mgDslashP QOP_D_mgDslashP
#define QOP_mgDslashPAll QOP_D_mgDslashPAll
#define QOP_mgDslashEo QOP_D_mgDslashEo
#define QOP_mgDslashEoProject QOP_D_mgDslashEoProject
#define QOP_mgDslashEoReconstruct QOP_D_mgDslashEoReconstruct
#define QOP_mgDslashEoReconstructP QOP_D_mgDslashEoReconstructP
#define QOP_Gcr QOP_D_Gcr
#define QOP_GcrSolveArgs QOP_D_GcrSolveArgs
#define QOP_gcrInit QOP_D_gcrInit
#define QOP_gcrFree QOP_D_gcrFree
#define QOP_gcrSet QOP_D_gcrSet
#define QOP_gcrSolve QOP_D_gcrSolve
#define QOP_gcrSolve2 QOP_D_gcrSolve2
#define QOP_gcrlsSolve QOP_D_gcrlsSolve
#define QOP_gcrcgSolve QOP_D_gcrcgSolve
#define QOP_gcrSolveA QOP_D_gcrSolveA
#define QOP_gcrSolveEo QOP_D_gcrSolveEo
#define QOP_gcrcgSolveA QOP_D_gcrcgSolveA
#define QOP_Cgls QOP_D_Cgls
#define QOP_CglsSolveArgs QOP_D_CglsSolveArgs
#define QOP_cglsInit QOP_D_cglsInit
#define QOP_cglsFree QOP_D_cglsFree
#define QOP_cglsSet QOP_D_cglsSet
#define QOP_cglsSolve QOP_D_cglsSolve
#define QOP_cgSolve QOP_D_cgSolve
#define QOP_cglsSolveA1 QOP_D_cglsSolveA1
#define QOP_cglsSolveA QOP_D_cglsSolveA
#define QOP_cglsSolveEo QOP_D_cglsSolveEo
#define QOP_cglsSolveEo1 QOP_D_cglsSolveEo1
#define QOP_cgSolveA QOP_D_cgSolveA
#define QOP_Bicgstab QOP_D_Bicgstab
#define QOP_BicgstabSolveArgs QOP_D_BicgstabSolveArgs
#define QOP_bicgstabInit QOP_D_bicgstabInit
#define QOP_bicgstabFree QOP_D_bicgstabFree
#define QOP_bicgstabSet QOP_D_bicgstabSet
#define QOP_bicgstabSolve QOP_D_bicgstabSolve
#define QOP_bicgstabSolveS QOP_D_bicgstabSolveS
#define QOP_bicgstabSolveA QOP_D_bicgstabSolveA
#define QOP_bicgstabSolveEo QOP_D_bicgstabSolveEo
#define QOP_MgVcycleArgs QOP_D_MgVcycleArgs
#define QOP_mgVcycle QOP_D_mgVcycle
#define QOP_eignhVN QOP_D_eignhVN
#define QOP_svdVN QOP_D_svdVN

#endif
