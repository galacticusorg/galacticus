/* Simple set of functions to return datatype IDs that are (annoyingly....) inaccessible from the HDF5 Fortran API. */

#include "H5Tpublic.h"

/* Return an ID for the basic C character string type */
hid_t H5T_C_S1_Get(void)
{
  return H5T_C_S1;
}

size_t H5T_Variable_Get(void)
{
  return H5T_VARIABLE;
}
