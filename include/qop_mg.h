typedef struct QOP_WilsonMgStruct QOP_WilsonMg;

QOP_WilsonMg *QOP_wilsonMgNew(void);
void QOP_wilsonMgFree(QOP_WilsonMg *wmg);
void QOP_wilsonMgSet(QOP_WilsonMg *wmg, int l, char *s, double val);
void QOP_wilsonMgSetArray(QOP_WilsonMg *wmg, int l, char *s, double *vals, int nval);
void QOP_wilsonMgSetLinks(QOP_WilsonMg *wmg, QOP_F3_FermionLinksWilson *wil);
void QOP_wilsonMgSetup(QOP_WilsonMg *wmg);

void QOP_F3_wilsonMgSolve(QOP_info_t *info, QOP_WilsonMg *wmg,
			  QOP_F3_FermionLinksWilson *flw,
			  QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg,
			  QLA_F_Real kappa, QDP_F3_DiracFermion *out,
			  QDP_F3_DiracFermion *in);
void QOP_D3_wilsonMgSolve(QOP_info_t *info, QOP_WilsonMg *wmg,
			  QOP_D3_FermionLinksWilson *flw,
			  QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg,
			  QLA_D_Real kappa, QDP_D3_DiracFermion *out,
			  QDP_D3_DiracFermion *in);

#if QOP_Precision == 'F'
#define QOP_wilsonMgSolve QOP_F3_wilsonMgSolve
#else
#define QOP_wilsonMgSolve QOP_D3_wilsonMgSolve
#endif
