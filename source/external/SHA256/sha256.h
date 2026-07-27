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

#ifndef SHA256_H
#define SHA256_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef struct{
    uint64_t size;        // Size of input in bytes
    uint32_t state[8];    // Current accumulation of hash
    uint8_t  input[64];   // Input to be used in the next step
    size_t   inputLength; // Number of bytes currently held in `input`
}SHA256Context;

void sha256Init    (SHA256Context *ctx                                     );
void sha256Update  (SHA256Context *ctx, const uint8_t *input, size_t length);
void sha256Finalize(SHA256Context *ctx,       uint8_t *digest              );

// Compute the SHA-256 digest of the file `fileName`, writing it to `hash` as a NULL-terminated string of 64 lowercase
// hexadecimal digits (so `hash` must have room for 65 characters). Returns 0 on success, or 1 if the file could not be read.
int sha256File(const char *fileName, char *hash);

#endif
