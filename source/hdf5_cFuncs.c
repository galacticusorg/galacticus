/* Simple set of interfaces to HDF5 C API functions. */

#include "H5Tpublic.h"

/* Call the H5close() function to completely shut down HDF5 */
void H5Close_C(void)
{
  herr_t error;  
  error=H5close();
  return;
}
