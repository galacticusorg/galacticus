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

// Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
//           2019
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

//% Implements Fortran-callable wrappers around the Linux mutex functions.

#include <stdlib.h>
#include <pthread.h>

int mutex_init(pthread_mutex_t **mutex, int recursive) {
  //% Fortran-callable wrapper around the pthread_mutex_init() function.
  int status;
  *mutex = malloc(sizeof(pthread_mutex_t));
  if (*mutex == NULL)
    return 1;
  if ( recursive ) {
    pthread_mutexattr_t attr;
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    status = pthread_mutex_init(*mutex, &attr);
    pthread_mutexattr_destroy(&attr);
  } else {
    status = pthread_mutex_init(*mutex, NULL);
  }
  return status;
}

int mutex_destroy(pthread_mutex_t *mutex) {
  //% Fortran-callable wrapper around the pthread_mutex_destroy() function.
  int status;
  status = pthread_mutex_destroy(mutex);
  if ( status != 0 )
    return status;
  free(mutex);
  return 0;
}

int* mutex_loc(pthread_mutex_t *mutex) {
  //% Return the address of a mutex.
  return (int*) mutex;
}
