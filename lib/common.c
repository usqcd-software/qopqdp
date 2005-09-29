#include <qop_internal.h>
#include <stdlib.h>
#include <sys/time.h>

double
dclock(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}
