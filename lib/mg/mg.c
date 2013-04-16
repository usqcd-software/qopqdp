#include <qop_internal.h>

#ifdef USE_MG

#if QOP_Precision == 'F'
#define QOPP(x) QOP_F_##x
#define QCDP(x) QOP_F_##x
#define QCDPC(x) QOP_F3_##x
#define QDPN(x) QDP_FN_##x
#else
#define QOPP(x) QOP_D_##x
#define QCDP(x) QOP_D_##x
#define QCDPC(x) QOP_D3_##x
#define QDPN(x) QDP_DN_##x
#endif

#define printf0 QOP_printf0

static int
hyper(QDP_Lattice *lat, int x[], void *args)
{
  int *bs = (int *)args;
  int nd = QDP_ndim_L(lat);
  int l[nd];
  QDP_latsize_L(lat, l);
  int r = 0;
  for(int i=nd-1; i>=0; i--) {
    r = r*(l[i]/bs[i]) + (x[i]/bs[i]);
  }
  return r;
}

#if 0

static void
f2c_map(int c[], int f[], QDP_Lattice *clat, QDP_Lattice *flat)
{
  int nd = QDP_ndim_L(clat);
  for(int i=0; i<nd; i++) {
    c[i] = (f[i]*QDP_coord_size_L(clat,i))/QDP_coord_size_L(flat,i);
  }
}

static void
c2f_map(int fmin[], int fmax[], int c[], QDP_Lattice *flat, QDP_Lattice *clat)
{
  int nd = QDP_ndim_L(flat);
  for(int i=0; i<nd; i++) {
    int fl = QDP_coord_size_L(flat,i);
    int cl = QDP_coord_size_L(clat,i);
    fmin[i] = (c[i]*fl+cl-1)/cl;
    fmax[i] = (c[i]*fl+fl-1)/cl;
  }
  if(1) {
#define PV(v,n) for(int l=0; l<n; l++) printf0(" %i", v[l]); printf0(" |")
#define CHKEQ for(int k=0; k<nd; k++) if(c[k]!=c2[k]) { printf0("error %s: %i %i %i\n", __func__, k, c[k], c2[k]); printf0("c,f0,f1 "); PV(c,nd); PV(fmin,nd); PV(fmax,nd); printf0("\n"); QDP_abort(1); }
    int c2[nd], f2[nd];
    f2c_map(c2, fmin, clat, flat);
    CHKEQ;
    f2c_map(c2, fmax, clat, flat);
    CHKEQ;
    for(int i=0; i<nd; i++) {
      for(int j=0; j<nd; j++) f2[j] = fmin[j];
      f2[i] = (f2[i]+QDP_coord_size_L(flat,i)-1)%QDP_coord_size_L(flat,i);
      f2c_map(c2, f2, clat, flat);
      c2[i] = (c2[i]+1)%QDP_coord_size_L(clat,i);
      CHKEQ;
    }
    for(int i=0; i<nd; i++) {
      for(int j=0; j<nd; j++) f2[j] = fmax[j];
      f2[i] = (f2[i]+1)%QDP_coord_size_L(flat,i);
      f2c_map(c2, f2, clat, flat);
      c2[i] = (c2[i]+QDP_coord_size_L(clat,i)-1)%QDP_coord_size_L(clat,i);
      CHKEQ;
    }
  }
}

#else

static void
f2c_map(int c[], int f[], int *lsc, int *lsf, int nd)
{
  for(int i=0; i<nd; i++) {
    c[i] = (f[i]*lsc[i])/lsf[i];
  }
}

static void
c2f_map(int fmin[], int fmax[], int c[], int *lsf, int *lsc, int nd)
{
  for(int i=0; i<nd; i++) {
    fmin[i] = (c[i]*lsf[i]+lsc[i]-1)/lsc[i];
    fmax[i] = (c[i]*lsf[i]+lsf[i]-1)/lsc[i];
  }
  if(1) {
#define PV(v,n) for(int l=0; l<n; l++) printf0(" %i", v[l]); printf0(" |")
#define CHKEQ for(int k=0; k<nd; k++) if(c[k]!=c2[k]) { printf0("error %s: %i %i %i\n", __func__, k, c[k], c2[k]); printf0("c,f0,f1 "); PV(c,nd); PV(fmin,nd); PV(fmax,nd); printf0("\n"); QDP_abort(1); }
    int c2[nd], f2[nd];
    f2c_map(c2, fmin, lsc, lsf, nd);
    CHKEQ;
    f2c_map(c2, fmax, lsc, lsf, nd);
    CHKEQ;
    for(int i=0; i<nd; i++) {
      for(int j=0; j<nd; j++) f2[j] = fmin[j];
      f2[i] = (f2[i]+lsf[i]-1)%lsf[i];
      f2c_map(c2, f2, lsc, lsf, nd);
      c2[i] = (c2[i]+1)%lsc[i];
      CHKEQ;
    }
    for(int i=0; i<nd; i++) {
      for(int j=0; j<nd; j++) f2[j] = fmax[j];
      f2[i] = (f2[i]+1)%lsf[i];
      f2c_map(c2, f2, lsc, lsf, nd);
      c2[i] = (c2[i]+lsc[i]-1)%lsc[i];
      CHKEQ;
    }
  }
}

#endif

static int
checkLocal(QDP_Lattice *flat, QDP_Lattice *clat)
{
  int nd = QDP_ndim_L(flat);
  int nsc = QDP_sites_on_node_L(clat);
  int xf0[nd], xf1[nd], xc[nd];
  int lsf[nd], lsc[nd];
  QDP_latsize_L(flat, lsf);
  QDP_latsize_L(clat, lsc);
  double nonlocal = 0;
  for(int i=0; i<nsc; i++) {
    QDP_get_coords_L(clat, xc, QDP_this_node, i);
    c2f_map(xf0, xf1, xc, lsf, lsc, nd);
    //printf0("xf0: "); PV(xf0,ndf); printf0("xf1: "); PV(xf1,ndf); printf0("\n");
    if(QDP_node_number_L(flat, xf0)!=QDP_this_node) { nonlocal = 1; break; }
    if(QDP_node_number_L(flat, xf1)!=QDP_this_node) { nonlocal = 1; break; }
  }
  //nonlocal = 1;
  QMP_sum_double(&nonlocal);
  int local = (nonlocal==0) ? 1 : 0;
  printf0("local = %i\n", local);
  return local;
}

static int
getpar(int *x, int nd)
{
  int p = 0;
  for(int i=0; i<nd; i++) p += x[i];
  return p&1;
}

QOP_MgBlock *
QOP_mgCreateBlockFromLattice(QDP_Lattice *fine, QDP_Lattice *coarse)
{
  QOP_MgBlock *mgb;
  QOP_malloc(mgb, QOP_MgBlock, 1);
  int nd = QDP_ndim_L(fine);
  int lsf[nd], lsc[nd];
  QDP_latsize_L(fine, lsf);
  QDP_latsize_L(coarse, lsc);

  mgb->fine = fine;
  mgb->coarse = coarse;

  mgb->local = 0;
  if(checkLocal(fine, coarse)) {
    mgb->local = 1;
    mgb->nlb = QDP_sites_on_node_L(coarse);
    QOP_malloc(mgb->lb, QOP_MgLocalBlock, mgb->nlb);
    int nsf = QDP_sites_on_node_L(fine);
    int *sites;
    QOP_malloc(sites, int, nsf);
    int xf[nd], xc[nd];
    for(int b=0; b<(mgb->nlb); b++) {
      for(int par=0; par<2; par++) {
	int n=0;
	mgb->lb[b].sites[par] = sites;
	for(int j=0; j<nsf; j++) {
	  //printf0("%i %i %i %p\n", QDP_this_node, j, QDP_sites_on_node_L(fine), xf);
	  QDP_get_coords_L(fine, xf, QDP_this_node, j);
	  //TRACE;
	  if(par==getpar(xf, nd)) {
	    f2c_map(xc, xf, lsc, lsf, nd);
	    if(QDP_index_L(coarse, xc)==b) {
	      sites[n] = j;
	      //printf0(" %i sites[%i] = %i\n", b, n, j);
	      n++;
	    }
	  }
	}
	mgb->lb[b].nsites[par] = n;
	sites += n;
      }
    }
  } else {
    // create subsets
    int ns = QDP_volume_L(coarse);
    int block[nd], ones[nd];
    for(int i=0; i<nd; i++) {
      block[i] = QDP_coord_size_L(fine, i)/QDP_coord_size_L(coarse, i);
      ones[i] = 1;
    }
    // need to fix this to use lattice size, not block size
    QDP_Subset *fs = QDP_create_subset_L(fine, hyper, block, nd*sizeof(int), ns);
    QDP_Subset *cs = QDP_create_subset_L(coarse, hyper, ones, nd*sizeof(int), ns);
    mgb->fs = fs;
    mgb->cs = cs;
    mgb->ns = ns;
  }

  return mgb;
}

QOP_MgBlock *
QOP_mgCreateBlock(QDP_Lattice *fine, int *block)
{
  // find coarse lattice dimensions
  int nd = QDP_ndim_L(fine);
  int fl[nd];
  QDP_latsize_L(fine, fl);
  int cl[nd];
  for(int i=0; i<nd; i++) {
    int k = fl[i]/block[i];
    if(k*block[i]!=fl[i]) {
      printf("error: fine lattice dim %i (%i) not divisible by block (%i)\n", i, fl[i], block[i]);
      QDP_abort(1);
    }
    cl[i] = k;
  }

  // create coarse lattice
  QDP_Lattice *coarse = QDP_create_lattice(QDP_get_default_layout(), NULL, nd, cl);
  return QOP_mgCreateBlockFromLattice(fine, coarse);
}

void
QOP_freeBlock(QOP_MgBlock *mgb)
{
  if(mgb->local) {
    QOP_free(mgb->lb[0].sites);
    QOP_free(mgb->lb);
  } else {
    QDP_destroy_subset(mgb->fs);
    QDP_destroy_subset(mgb->cs);
  }
  QDP_destroy_lattice(mgb->coarse);
  QOP_free(mgb);
}

#endif
