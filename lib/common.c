#include <stdlib.h>
#include <sys/time.h>
#include <qop_internal.h>

QOP_common_t QOP_common = {.inited=0};


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
    QDP_set_latsize(layout->latdim, layout->latsize);
    QDP_create_layout();
    QOP_common.we_inited_qdp = 1;
  }
  if( (layout->machdim>=0) && (QMP_logical_topology_is_declared()) ) {
    int i, n, qmpndim, error=0;
    const int *qmpdims;
    n = layout->machdim;
    qmpndim = QMP_get_logical_number_of_dimensions();
    if(qmpndim>n) n = qmpndim;
    qmpdims = QMP_get_logical_dimensions();
    for(i=0; i<n; i++) {
      int qopsize=1, qmpsize=1;
      if(i<layout->machdim) qopsize = layout->machsize[i];
      if(i<qmpndim) qmpsize = qmpdims[i];
      if(qopsize!=qmpsize) error = 1;
    }
    if(error) {
      if(QDP_this_node==0) {
	printf("QOP Warning: QOP machine != QMP logical topology\n");
	printf("QOP machsize =");
	for(i=0; i<layout->machdim; i++) printf(" %i", layout->machsize[i]);
	printf("\nQMP topology =");
	for(i=0; i<qmpndim; i++) printf(" %i", qmpdims[i]);
	printf("\n");
	fflush(stdout);
      }
    }
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
