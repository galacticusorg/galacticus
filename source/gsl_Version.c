/* Function that returns the GSL version string. */

#include "gsl/gsl_version.h"

/* Return the GSL version string. */
const void * GSL_Get_Version(void)
{
  return gsl_version;
}
