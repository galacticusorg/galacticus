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

//% Implements Fortran-callable wrappers around the GNU C Library regular expression functions.

#include <regex.h>
#include <stdlib.h>
#include <stdio.h>

regex_t * Regular_Expression_Construct_C(const char *pattern) {
  //% Fortran-callable wrapper around the GNU C Library {\tt regcomp} regular expression compiling function.
  regex_t *r;
  int s;
  r=(regex_t *)malloc(sizeof(*r));
  s=regcomp(r,pattern,REG_EXTENDED);
  if (s != 0) {
    printf("Regular_Expression_Construct_C() failed to compile [%s]\n",pattern);
    abort();
  }
  return r;
}

void Regular_Expression_Destruct_C(regex_t *r) {
  //% Fortran-callable function to free a {\tt regex_t} struct.
  regfree(r);
  free(r);
}

int Regular_Expression_Match_C(regex_t *r, const char *string) {
  //% Fortran-callable wrapper around the GNU C Library {\tt regexec} regular expression matching function.
 return regexec (r,string,0,NULL,0);
}
