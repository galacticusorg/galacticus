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

//% Implements Fortran-callable wrappers around errno, for reporting system errors.

#include <stdio.h>
#include <errno.h>
#include <string.h>

int errnoGet_C() {
  //% Return the value of errno for the calling thread.
  return errno;
}

void errnoClear_C() {
  //% Set errno to zero for the calling thread.
  errno = 0;
}

void errnoDescription_C(int errorNumber, char *description, int descriptionLength) {
  //% Write a null-terminated description of the given error number into the supplied buffer. Uses strerror_r(), as strerror() is not thread safe.
  if (descriptionLength <= 0) return;
  description[0] = '\0';
#if defined(_GNU_SOURCE) && defined(__GLIBC__)
  /* GNU strerror_r() returns a pointer to the description, which may or may not be the supplied buffer. */
  const char *descriptionGNU = strerror_r(errorNumber, description, (size_t)descriptionLength);
  if (descriptionGNU != description) snprintf(description, (size_t)descriptionLength, "%s", descriptionGNU);
#else
  /* XSI strerror_r() writes into the supplied buffer, and returns non-zero on failure. */
  if (strerror_r(errorNumber, description, (size_t)descriptionLength) != 0) snprintf(description, (size_t)descriptionLength, "unknown error %d", errorNumber);
#endif
}
