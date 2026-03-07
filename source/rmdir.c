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

//% Implements Fortran-callable wrappers around the Linux rmdir function.

#include <unistd.h>
#include <errno.h>
#include <stdio.h>

int rmdir_C(const char *path) {
  //% Fortran-callable wrapper around the rmdir() function to remove a directory.
  int status;
  status = rmdir(path);
  if ( status == -1 ) {
    int err = errno;
    if ( err == ENOENT ) {
      /* Directory does not exist - this is acceptable */
      return 0;
    } else {
      if        ( err == EACCES       ) {
	printf("rmdir failed: write access to the directory containing path was not allowed\n");
      } else if ( err == EBUSY        ) {
	printf("rmdir failed: path is currently in use by the system or some process that prevents its removal\n");
      } else if ( err == EFAULT       ) {
	printf("rmdir failed: path points outside your accessible address space\n");
      } else if ( err == EINVAL       ) {
	printf("rmdir failed: path has .  as last component\n");
      } else if ( err == ELOOP        ) {
	printf("rmdir failed: too many symbolic links were encountered in resolving path\n");
      } else if ( err == ENAMETOOLONG ) {
	printf("rmdir failed: path was too long\n");
      } else if ( err == ENOENT       ) {
	printf("rmdir failed: a directory component in path does not exist or is adangling symbolic link\n");
      } else if ( err == ENOMEM       ) {
	printf("rmdir failed: insufficient kernel memory was available\n");
      } else if ( err == ENOTDIR      ) {
	printf("rmdir failed: path, or a component used as a directory in path, is not, in fact, a directory\n");
      } else if ( err == ENOTEMPTY    ) {
	printf("rmdir failed: path contains entries other than .  and ..; or, path has .. as its final component\n");
      } else if ( err == EPERM        ) {
	printf("rmdir failed: permissions issue\n");
      } else if ( err == EROFS        ) {
	printf("rmdir failed: path refers to a directory on a read-only filesystem\n");
      }
      return err;
    }
  } else {
    return status;
  }
}
