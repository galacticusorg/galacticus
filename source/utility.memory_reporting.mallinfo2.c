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

//% Implements Fortran-callable wrappers around the Linux mallinfo2() function.

#ifdef __APPLE__
#include <stdlib.h>
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#include <unistd.h>

#ifdef __linux__
#ifdef __GNUC__
#include <gnu/libc-version.h>
#ifdef __GLIBC__
#if ( __GLIBC__ == 2 && __GLIBC_MINOR__ >= 33 ) || __GLIBC__ > 2
#define mallinfo2_available
#endif
#endif
#endif
#endif

size_t mallinfo2_c() {
  //% Fortran-callable wrapper around the mallinfo2() function to get memory usage information.
#ifdef __APPLE__
  struct mstats ms = mstats();
  size_t uordblks = ms.bytes_used;
  return uordblks;
#else
#ifdef mallinfo2_available
  struct mallinfo2 info;
  info = mallinfo2();
  return info.uordblks+info.usmblks+info.hblkhd;
#else
  struct mallinfo info;
  info = mallinfo();
  size_t uordblks = info.uordblks;
  size_t usmblks  = info.usmblks ;
  size_t hblkhd   = info.hblkhd  ;
  size_t one      = 1            ;
  // Old mallinfo() uses ints which can overflow and become negative. Trap that here - there's nothing more we can do about
  // this. Newer glibcs support mallinfo2() which does not have this problem.
  if (uordblks < 0 || uordblks > one<<32) uordblks=0;
  if (usmblks  < 0 || usmblks  > one<<32) usmblks =0;
  if (hblkhd   < 0 || hblkhd   > one<<32) hblkhd  =0;
  return uordblks+usmblks+hblkhd;
#endif
#endif
}

size_t gettotalsystemmemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
