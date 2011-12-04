#include <stdlib.h>
#include <sys/time.h>
#include <qop_internal.h>

QOP_common_t QOP_common = {.inited=0};
extern QDP_Layout *QOP_layout_user;

static int
compare_sizes(int n1, int *v1, char *s1, int n2, int *v2, char *s2)
{
  int i, n, error=0;

  n = n1;
  if(n2>n1) n = n2;
  for(i=0; i<n; i++) {
    int t1=1, t2=1;
    if(i<n1) t1 = v1[i];
    if(i<n2) t2 = v2[i];
    if(t1!=t2) error = 1;
  }
  if(error) {
    if(QDP_this_node==0) {
      printf("QOP Warning: %s != %s\n", s1, s2);
      printf("%s =", s1);
      for(i=0; i<n1; i++) printf(" %i", v1[i]);
      printf("\n%s =", s2);
      for(i=0; i<n2; i++) printf(" %i", v2[i]);
      printf("\n");
      fflush(stdout);
    }
  }
  return error;
}

  /*********************/
  /*  Public routines  */
  /*********************/

QOP_status_t
QOP_init(QOP_layout_t *layout)
{
  QOP_status_t retval = QOP_SUCCESS;

  CHECK_NOT_INIT;
  QOP_common.inited = 1;
  QOP_common.verbosity = 0;
  QOP_common.proflevel = 0;
  QOP_common.we_inited_qdp = 0;
  QOP_common.ndim = layout->latdim;

  if(!QDP_is_initialized()) {
    QDP_initialize(NULL, NULL);
    QDP_set_default_layout(QOP_layout_user);
    QDP_set_latsize(layout->latdim, layout->latsize);
    QDP_create_layout();
    QOP_common.we_inited_qdp = 1;
  } else {
    int error;
    int qdplatdim = QDP_ndim();
    int qdplatsize[qdplatdim];
    QDP_latsize(qdplatsize);
    error = compare_sizes(layout->latdim, layout->latsize, "QOP lattice",
			  qdplatdim, qdplatsize, "QDP lattice");
    if(error) retval = QOP_FAIL;
  }

  if( (layout->machdim>=0) && (QMP_logical_topology_is_declared()) ) {
    int error;
    int qmpndim;
    const int *qmpdims;
    qmpndim = QMP_get_logical_number_of_dimensions();
    qmpdims = QMP_get_logical_dimensions();
    error = compare_sizes(layout->machdim, layout->machsize, "QOP machsize",
			  qmpndim, (int *)qmpdims, "QMP topology");
    if(error) retval = QOP_FAIL;
  }

  return retval;
}

QOP_status_t
QOP_finalize(void)
{
  CHECK_INIT;
  QOP_common.inited = 0;
  if(QOP_common.we_inited_qdp) {
    QDP_finalize();
  }
  return QOP_SUCCESS;
}

int
QOP_is_initialized(void)
{
  return QOP_common.inited;
}

int
QOP_verbose(int level)
{
  int old = QOP_common.verbosity;
  QOP_common.verbosity = level;
  return old;
}

int
QOP_profcontrol(int level)
{
  int old = QOP_common.proflevel;
  QOP_common.proflevel = level;
  return old;
}

int
QOP_node_number_raw(int coords[])
{
  CHECK_INIT;
  return QDP_node_number(coords);
}

int
QOP_node_index_raw_V(int coords[], QOP_evenodd_t evenodd)
{
  CHECK_INIT;
  return QDP_index(coords);
}

int
QOP_node_index_raw_D(int coords[], QOP_evenodd_t evenodd)
{
  CHECK_INIT;
  return QDP_index(coords);
}

int
QOP_node_index_raw_G(int coords[], QOP_evenodd_t evenodd)
{
  CHECK_INIT;
  return QDP_index(coords);
}

int
QOP_node_index_raw_F(int coords[], QOP_evenodd_t evenodd)
{
  CHECK_INIT;
  return QDP_index(coords);
}

int
QOP_sites_on_node_raw_V(QOP_evenodd_t evenodd)
{
  return QDP_sites_on_node;
}

int
QOP_sites_on_node_raw_D(QOP_evenodd_t evenodd)
{
  return QDP_sites_on_node;
}

int
QOP_sites_on_node_raw_G(QOP_evenodd_t evenodd)
{
  return QDP_sites_on_node;
}

int
QOP_sites_on_node_raw_F(QOP_evenodd_t evenodd)
{
  return QDP_sites_on_node;
}


  /***********************/
  /*  Internal routines  */
  /***********************/

double
QOP_time(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}
