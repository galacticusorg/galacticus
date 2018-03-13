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

//% Implements Fortran-callable wrappers around the Linux file locking functions.

#include <flock_config.h>
#ifdef OFDAVAIL
#define _GNU_SOURCE
#endif
#include <unistd.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

struct lockDescriptor {
  struct flock fl;
  int fd;
  char *name;
};

void flock_C(const char *name, struct lockDescriptor **ld, int lockIsShared) {
  //% Fortran-callable wrapper around the fcntl() function to lock a file.
  int e;
  *ld = malloc(sizeof(struct lockDescriptor));
  size_t destination_size = strlen(name);
  (*ld)->name = malloc(sizeof(char)*(destination_size+1));
  (*ld)->name[destination_size] = '\0';
  strncpy((*ld)->name,name,destination_size);
  if (lockIsShared == 0) {
    (*ld)->fl.l_type  =F_WRLCK;
  } else {
    (*ld)->fl.l_type  =F_RDLCK;
  }
  (*ld)->fl.l_whence=SEEK_SET;
  (*ld)->fl.l_start =0;        /* Offset from l_whence                        */
  (*ld)->fl.l_len   =1;        /* length, 0 = to EOF                          */
#ifdef OFDAVAIL
  (*ld)->fl.l_pid   =0;        /* no PID for Linux open file descriptor locks */
#else
  (*ld)->fl.l_pid   =getpid(); /* our PID                                     */
#endif
  (*ld)->fd=open(name,O_RDWR | O_CREAT,S_IRUSR | S_IWUSR);
  if ( (*ld)->fd < 0 ) {
    printf("flock_C(): opening file '%s' failed\n",name);
    printf("  -> error no is : %d\n", errno);
    printf("  -> error description is : %s\n",strerror(errno));
    abort();
  }
#ifdef OFDAVAIL
  if (fcntl((*ld)->fd, F_OFD_SETLKW, &(*ld)->fl) == -1) {
#else
  if (fcntl((*ld)->fd, F_SETLKW    , &(*ld)->fl) == -1) {
#endif
    if (errno == ENOSYS) {
      /* File locking is not implemented on this system. Fail silently, returning a suitable error code in the lock descriptor */
      (*ld)->fd = -1;
    } else {    
      if (errno == EBADF) {
	printf("flock_C(): bad file descriptor [EBADF]: %s\n",name);
      } else if (errno == EINTR) {
	printf("flock_C(): [EINTR]\n");
      } else if (errno == EINVAL) {
      printf("flock_C(): [EINVAL]\n");
      } else if (errno == ENOLCK) {
	printf("flock_C(): [ENOLCK]\n");
      } else if (errno == EOVERFLOW) {
	printf("flock_C(): [EOVERFLOW]\n");
      } else if (errno == EDEADLK) {
	printf("flock_C(): [EDEADLK]\n");
      } else {
	printf("flock_C(): unknown error [%d]\n",errno);
      }
      abort();
    }
  }
}

void funlock_C(struct lockDescriptor **ld) {
  //% Fortran-callable wrapper around the fcntl() function to unlock a file.
  /* Skip unlock if the original lock failed due to locking not being available on the system */
  if ( (*ld)->fd != -1 ) {
    (*ld)->fl.l_type  =F_UNLCK;
    close((*ld)->fd);
#ifdef OFDAVAIL
    fcntl((*ld)->fd, F_OFD_SETLK, &(*ld)->fl); /* set the region to unlocked */
#else
    fcntl((*ld)->fd, F_SETLK    , &(*ld)->fl); /* set the region to unlocked */
#endif
    free((*ld)->name);
    free(*ld);
  }
  return;
}
