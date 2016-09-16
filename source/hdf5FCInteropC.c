/* C-side determination of interoperable HDF5 types  */

#include "H5Tpublic.h"
#include <stdio.h>

int main() {
  /* Determine the size of the "herr_t" type */
  if        (sizeof(herr_t)==sizeof(int      )) {
    printf(" herr_t = c_int\n"        );
  } else if (sizeof(herr_t)==sizeof(short    )) {
    printf(" herr_t = c_short\n"      );    
  } else if (sizeof(herr_t)==sizeof(long     )) {
    printf(" herr_t = c_long\n"       );    
  } else if (sizeof(herr_t)==sizeof(long long)) {
    printf(" herr_t = c_long_long\n");    
  } else if (sizeof(herr_t)==sizeof(size_t   )) {
    printf(" herr_t = c_size_t\n"     );
  } else {
    printf(" herr_t = ?????\n"        );
  }
}
