/* Simple set of interfaces to HDF5 C API functions. */

#include "H5Tpublic.h"

/* Return an ID for the basic C character string type */
void H5Close_C(void)
{
  herr_t error;  
  error=H5close();
  return;
}
