#include <stdlib.h>
#include <stdio.h>
#include <qop_internal.h>

// prevent use of default lattice variables/functions
#define QDP_sites_on_node ERROR
#define QDP_ndim() ERROR
#define QDP_coord_size(i) ERROR
#define QDP_latsize(ls) ERROR
#define QDP_volume() ERROR
#define QDP_numsites(node) ERROR
#define QDP_node_number(x) ERROR
#define QDP_index(x) ERROR
#define QDP_get_coords(x, node, index) ERROR
#define QDP_all ERROR
#define QDP_even_and_odd ERROR
#define QDP_even ERROR
#define QDP_odd ERROR
#define QDP_neighbor ERROR
////////////////////////////////////

typedef struct {
  int *len;          /* lattice dimensions */
  int *nsquares;     /* number of hypercubes in each direction */
  int ndim;
  int numsites;
} params;

static int prime[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
#define MAXPRIMES (sizeof(prime)/sizeof(int))

static void
get_lex_x(int *x, int l, int *s, int ndim)
{
  //for(int i=0; i<ndim; i++) {
  for(int i=ndim-1; i>=0; --i) {
    x[i] = l % s[i];
    l = l / s[i];
  }
}

#if 0
static int
get_lex_i(int *x, int *s, int ndim)
{
  int l = 0;
  //for(int i=ndim-1; i>=0; --i) {
  for(int i=0; i<ndim; i++) {
    if(x[i]>=s[i]) { l = -1; break; }
    l = l*s[i] + (x[i]%s[i]);
  }
  return l;
}
#endif

#if 0
static void
node2coord(int *x, int n, params *p)
{
  for(int i=0; i<p->ndim; i++) {
    x[i] = n % p->nsquares[i];
    n = n / p->nsquares[i];
  }
}

static int
coord2node(int *x, params *p)
{
  int l = 0;
  for(int i=p->ndim-1; i>=0; --i) {
    if(x[i]>=p->nsquares[i]) { l = -1; break; }
    l = l*p->nsquares[i] + (x[i]%p->nsquares[i]);
  }
  return l;
}
#endif

static void
layout_hyper_eo_setup(QDP_Lattice *lat, void *args)
{
  QDP_allocate_lattice_params(lat, sizeof(params));
  params *p = (params *) QDP_get_lattice_params(lat);

  int nd = QDP_ndim_L(lat);

  p->ndim = nd;
  p->len = (int *) malloc(nd*sizeof(int));
  p->nsquares = (int *) malloc(nd*sizeof(int));
  int *len = p->len;
  int *nsquares = p->nsquares;
  QDP_latsize_L(lat, len);

  if(QMP_logical_topology_is_declared()) { // use logical topology
    int nd2 = QMP_get_logical_number_of_dimensions();
    const int *nsquares2 = QMP_get_logical_dimensions();
    for(int i=0; i<nd; i++) {
      if(i<nd2) nsquares[i] = nsquares2[i];
      else nsquares[i] = 1;
    }
  } else if(QMP_get_msg_passing_type()!=QMP_SWITCH) {
    int nd2 = QMP_get_allocated_number_of_dimensions();
    const int *nsquares2 = QMP_get_allocated_dimensions();
    for(int i=0; i<nd; i++) {
      if(i<nd2) nsquares[i] = nsquares2[i];
      else nsquares[i] = 1;
    }
  } else { /* not QMP_GRID */
    int *squaresize = (int *) malloc(nd*sizeof(int));
    int *extrafactors = (int *) malloc(nd*sizeof(int));
    for(int i=0; i<nd; ++i) {
      squaresize[i] = len[i];
      extrafactors[i] = 1;
      nsquares[i] = 1;
    }

    /* Figure out dimensions of rectangle */
    int n = QMP_get_number_of_nodes();   /* remaining number of nodes to be factored */
    int k = MAXPRIMES-1;
    while(n>1) {
      /* figure out which prime to divide by starting with largest */
      /* if no factor found, assume n is prime */
      while( (k>=0) && (n%prime[k]!=0) ) --k;
      int pfac = (k>=0) ? prime[k] : n;

      /* figure out which direction to divide */
      /* find largest divisible dimension of h-cubes */
      /* if one direction with largest dimension has already been
	 divided, divide it again.  Otherwise divide first direction
	 with largest dimension. */
      int j = -1;
      for(int i=0; i<nd; i++) {
	if(squaresize[i]%pfac==0) {
	  if( (j<0) || (extrafactors[j]*squaresize[i]>extrafactors[i]*squaresize[j]) ) {
	    j = i;
	  } else if(extrafactors[j]*squaresize[i]==extrafactors[i]*squaresize[j]) {
	    if((nsquares[j]==1)||(nsquares[i]!=1)) j = i;
	  }
	}
      }

      /* This can fail if we run out of prime factors in the dimensions */
      /* then just choose largest dimension */
      if(j<0) {
	for(int i=0; i<nd; i++) {
	  if( (j<0) || (extrafactors[j]*squaresize[i]>extrafactors[i]*squaresize[j]) ) {
	    j = i;
	  } else if(extrafactors[j]*squaresize[i]==extrafactors[i]*squaresize[j]) {
	    if((nsquares[j]==1)||(nsquares[i]!=1)) j = i;
	  }
	}
	n /= pfac;
	extrafactors[j] *= pfac;
	nsquares[j] *= pfac;
      } else {
	n /= pfac;
	squaresize[j] /= pfac;
	nsquares[j] *= pfac;
      }
    }
    free(squaresize);
    free(extrafactors);
  } /* not QMP_GRID */

  /* setup QMP logical topology */
#if 1
  if(!QMP_logical_topology_is_declared()) {
    QMP_declare_logical_topology(nsquares, nd);
  }
#endif

  int numsites = 1;
  int mcoord[nd];
  QMP_get_logical_coordinates_from2(mcoord, QDP_this_node);
  //node2coord(mcoord, QDP_this_node, p);
  for(int i=0; i<nd; ++i) {
    int x0 = (mcoord[i]*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
    int x1 = ((mcoord[i]+1)*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
    numsites *= x1-x0;
  }
  p->numsites = numsites;

  if(QDP_this_node==0) {
    printf("ndim = %i\n", p->ndim);
    printf("numsites = %i\n", p->numsites);
    printf("len =");
    for(int i=0; i<p->ndim; i++) printf(" %i", p->len[i]);
    printf("\n");
    printf("nsquares =");
    for(int i=0; i<p->ndim; i++) printf(" %i", p->nsquares[i]);
    printf("\n");
  }
}

static void
layout_hyper_eo_free(QDP_Lattice *lat)
{
  params *p = (params *) QDP_get_lattice_params(lat);
  free(p->len);
  free(p->nsquares);
}


static int
layout_hyper_eo_numsites(QDP_Lattice *lat, int node)
{
  params *p = (params *) QDP_get_lattice_params(lat);
  if(node==QDP_this_node) {
    return p->numsites;
  } else {
    int numsites = 1;
    int nd = p->ndim;
    int mcoord[nd];
    QMP_get_logical_coordinates_from2(mcoord, node);
    //node2coord(mcoord, node, p);
    for(int i=0; i<nd; ++i) {
      int x0 = (mcoord[i]*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
      int x1 = ((mcoord[i]+1)*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
      numsites *= x1-x0;
    }
    return numsites;
  }
}

static int
layout_hyper_eo_node_number(QDP_Lattice *lat, const int x[])
{
  params *p = (params *) QDP_get_lattice_params(lat);
  int i, m[p->ndim];

  for(i=0; i<p->ndim; i++) {
    m[i] = (x[i]*p->nsquares[i])/p->len[i];
  }
  return QMP_get_node_number_from(m);
  //return coord2node(m, p);
}

static int
layout_hyper_eo_index(QDP_Lattice *lat, const int x[])
{
  params *p = (params *) QDP_get_lattice_params(lat);
  int s=0, l=0;

  //for(int i=p->ndim-1; i>=0; --i) {
  for(int i=0; i<p->ndim; ++i) {
    int m = (x[i]*p->nsquares[i])/p->len[i];
    int x0 = (m*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
    int x1 = ((m+1)*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
    l = l*(x1-x0) + x[i]-x0;
    s += x[i];
  }

  if( s%2==0 ) { /* even site */
    l /= 2;
  } else {
    l = (l+p->numsites)/2;
  }
  return l;
}

static void
layout_hyper_eo_get_coords(QDP_Lattice *lat, int x[], int node, int index)
{
  params *p = (params *) QDP_get_lattice_params(lat);
  int nd = p->ndim;
  int m[nd], dx[nd], sx[nd];

  QMP_get_logical_coordinates_from2(m, node);
  //node2coord(m, node, p);

  int s0 = 0;
  int n0 = 1;
  for(int i=0; i<nd; ++i) {
    x[i] = (m[i]*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
    int x1 = ((m[i]+1)*p->len[i]+p->nsquares[i]-1)/p->nsquares[i];
    dx[i] = x1 - x[i];
    s0 += x[i];
    n0 *= dx[i];
  }

  int neven = (n0 + 1 - (s0&1))/2;
  //if(QDP_this_node==0) printf("neven = %i\n", neven);
  if(index<neven) {
    int l = 2*index;
    get_lex_x(sx, l, dx, nd);
    int s1 = s0;
    for(int i=0; i<nd; ++i) s1 += sx[i];
    if((s1&1)!=0) {
      get_lex_x(sx, l+1, dx, nd);
    }
  } else {
    int l = 2*index - n0 + ((n0&1)*(s0&1));
    get_lex_x(sx, l, dx, nd);
    int s1 = s0;
    for(int i=0; i<nd; ++i) s1 += sx[i];
    //if(QDP_this_node==0) printf("index = %i  n0 = %i  s0 = %i  l = %i  s1 = %i\n",index,n0,s0,l,s1);
    if((s1&1)==0) {
      get_lex_x(sx, l+1, dx, nd);
    }
  }
  for(int i=0; i<nd; ++i) x[i] += sx[i];

  if(QDP_index_L(lat, x)!=index) {
    if(QDP_this_node==0) {
      fprintf(stderr,"QDP: error in layout!\n");
      fprintf(stderr,"%i %i  -> ", node, index);
      for(int i=0; i<nd; i++) fprintf(stderr," %i", x[i]);
      fprintf(stderr,"  ->  %i %i\n", QDP_node_number_L(lat,x), QDP_index_L(lat,x));
    }
    QDP_abort(1);
    exit(1);
  }
}

static QDP_Layout layout_hyper_eo = {
  layout_hyper_eo_setup,
  layout_hyper_eo_free,
  layout_hyper_eo_numsites,
  layout_hyper_eo_node_number,
  layout_hyper_eo_index,
  layout_hyper_eo_get_coords
};
QDP_Layout *QOP_layout_user = &layout_hyper_eo;
