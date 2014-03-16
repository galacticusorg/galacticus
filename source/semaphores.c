// Copyright 2009, 2010, 2011 Andrew Benson <abenson@obs.carnegiescience.edu>
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

//% Implements Fortran-callable wrappers around the Linux semaphore functions.

#include <sys/stat.h>
#include <fcntl.h>
#include <semaphore.h>
#include <errno.h>
#include <stdlib.h>

sem_t * Semaphore_Open_C(const char *name, int initialValue) {
  //% Fortran-callable wrapper around the Linux semaphore open function.
  sem_t *s;
  s = sem_open(name,O_CREAT,S_IRUSR | S_IWUSR,initialValue);
  if ( s == SEM_FAILED ) {
    if ( errno == EACCES       ) printf("sem_open() failed: access\n"         );
    if ( errno == EEXIST       ) printf("sem_open() failed: exists\n"         );
    if ( errno == EINVAL       ) printf("sem_open() failed: incorrect value\n");
    if ( errno == ENAMETOOLONG ) printf("sem_open() failed: name too long\n"  );
    if ( errno == ENOENT       ) printf("sem_open() failed: no such path\n"   );
    if ( errno == ENOSPC       ) printf("sem_open() failed: no space\n"       );
    abort();
  }
  return s;
}

void Semaphore_Close_C(sem_t *s) {
  //% Fortran-callable wrapper around the Linux semaphore close function.
  int e;
  e = sem_close(s);
  if ( e == -1 ) {
    if ( errno == EINVAL ) printf("sem_close() failed: incorrect value\n");
    abort();
  }
}

void Semaphore_Wait_C(sem_t *s) {
  //% Fortran-callable wrapper around the Linux semaphore wait function.
  int e;
  e = -1;
  while ( e == -1 ) {
    e = sem_wait(s);
    if ( e == -1 ) {
      if ( errno == EINVAL ) {
	printf("sem_wait() failed: incorrect value\n");
	abort();
      }
      /* if ( errno == EINTR  ) printf("sem_wait() failed: interrupted\n"    ); */
    }
  }
}

void Semaphore_Post_C(sem_t *s) {
  //% Fortran-callable wrapper around the Linux semaphore post function.
  int e;
  e = sem_post(s);
  if ( e == -1 ) {
    if ( errno == EINVAL ) printf("sem_post() failed: incorrect value\n");
    abort();
  }
}

void Semaphore_Unlink_C(const char *name) {
  //% Fortran-callable wrapper around the Linux semaphore unlink function.
  int e;
  e = sem_unlink(name);
  if ( e == -1 ) {
    if ( errno == EACCES       ) printf("sem_unlink() failed: access\n"           );
    if ( errno == ENAMETOOLONG ) printf("sem_unlink() failed: name too long\n"    );
    if ( errno == ENOENT       ) printf("sem_unlink() failed: no such path\n"     );
    if ( errno == EFAULT       ) printf("sem_unlink() failed: incorrect address\n");
    abort();
  }
}
