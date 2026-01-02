// Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
//           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
//    Andrew Benson <abenson@carnegiescience.edu>
//
// This file is part of Galacticus.
//
//    Galacticus is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Galacticus is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

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
