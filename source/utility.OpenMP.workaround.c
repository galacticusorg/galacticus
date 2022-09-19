/* From glibc version 2.34 onward libpthread.so no longer exists and is part of the main libc.so.       */
/* This causes problems for statitically-linked OpenMP codes because (for some reason) the function     */
/* pthread_mutex_destroy() is not linked so results in a segfault when called. The following workaround */
/* adds a call to pthread_mutex_destroy() so that is is included in the static executable.              */

#ifdef STATIC
#ifdef __GNUC__
#include <gnu/libc-version.h>
#ifdef __GLIBC__
#if ( __GLIBC__ == 2 && __GLIBC_MINOR__ >= 34 ) || __GLIBC__ > 2
#include "pthread.h"
#define nullptr ((void*)0)

void pthread_workaround() {
  pthread_mutex_destroy((pthread_mutex_t *) nullptr);
}

#endif
#endif
#endif
#endif
