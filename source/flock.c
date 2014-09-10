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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/file.h>

int flock_C(const char *name) {
  //% Fortran-callable wrapper around the flock() function to lock a file.
  int fd;
  fd=open(name,O_RDONLY);
  return flock(fd,LOCK_EX);
}

void funlock_C(int fd) {
  //% Fortran-callable wrapper around the flock() function to unlock a file.
  flock(fd,LOCK_UN);
  close(fd);
  return;
}
