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

/* Simple set of functions to return datatype IDs that are (annoyingly....) inaccessible from the HDF5 Fortran API. */

#include <stddef.h>
#include "H5Tpublic.h"

/* Return an ID for the basic C character string type */
hid_t H5T_C_S1_Get(void)
{
  return H5T_C_S1;
}

size_t H5T_Variable_Get(void)
{
  return H5T_VARIABLE;
}

/* Length of the string fields inside unitType - must match unitStringLength in utility.units.F90 (Units_MetaData module) */
#define UNIT_STRING_LENGTH 512

/* C mirror of the Fortran unitType bind(C) derived type.
   The field order and sizes must match exactly. */
typedef struct {
  double unitsInSI;
  char   description[UNIT_STRING_LENGTH];
  char   quantity   [UNIT_STRING_LENGTH];
  int    isComoving;
} unitType_c;

/* Create and return an HDF5 compound datatype that matches unitType_c.
   The caller is responsible for closing the returned type ID with H5Tclose(). */
hid_t HDF5_Unit_Type_Create(void)
{
  hid_t  compound_id, string_id;
  herr_t err;

  /* Build a fixed-length C string type of UNIT_STRING_LENGTH bytes. */
  string_id = H5Tcopy(H5T_C_S1);
  err = H5Tset_size(string_id, (size_t)UNIT_STRING_LENGTH);
  err = H5Tset_strpad(string_id, H5T_STR_NULLTERM);
  err = H5Tset_cset(string_id, H5T_CSET_ASCII);

  /* Build the compound type. */
  compound_id = H5Tcreate(H5T_COMPOUND, sizeof(unitType_c));
  err = H5Tinsert(compound_id, "unitsInSI"  , offsetof(unitType_c, unitsInSI  ), H5T_NATIVE_DOUBLE);
  err = H5Tinsert(compound_id, "description", offsetof(unitType_c, description), string_id        );
  err = H5Tinsert(compound_id, "quantity"   , offsetof(unitType_c, quantity   ), string_id        );
  err = H5Tinsert(compound_id, "isComoving" , offsetof(unitType_c, isComoving ), H5T_NATIVE_INT   );

  H5Tclose(string_id);
  (void)err; /* suppress unused-variable warning */
  return compound_id;
}
