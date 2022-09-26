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

//% Implements Fortran-callable wrappers around the POSIX access() function.

#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>

int access_C(const char *name) {
  //% Fortran-callable wrapper around the POSIX access() function to check file existance.
  return access(name,F_OK);
}

void syncdir_C(const char *name) {
  //% Fortran-callable function to sync a directory. This is done by opening and closing the directory. According to the this
  //% \href{https://stackoverflow.com/a/30630912}{answer} on Stackoverflow this will invalidate the NFS cache for this directory.

  printf("symcdir_C() #1 '%s'\n",name);

  DIR *fd = opendir (name);

  printf("symcdir_C() #2 '%s'\n",name);

  int i   = closedir(  fd);
  printf("symcdir_C() #3 '%s'\n",name);
}
