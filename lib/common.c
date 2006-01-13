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
  CHECK_NOT_INIT;
  QOP_common.inited = 1;
  QOP_common.verbosity = 0;
  QOP_common.proflevel = 0;
  QOP_common.we_inited_qdp = 0;
  QOP_common.ndim = layout->latdim;
  //if(!QDP_initialized()) {
  if(0) {
    QDP_initialize(NULL, NULL);
    QDP_set_latsize(layout->latdim, layout->latsize);
    QDP_create_layout();
    QOP_common.we_inited_qdp = 1;
  }
  return QOP_SUCCESS;
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
