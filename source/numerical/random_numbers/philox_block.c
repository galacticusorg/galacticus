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

//+    Contributions to this file made by: Andrew Benson, Claude.

/*
 * The Philox4x64-10 counter-based random number generator of Salmon et al. (2011; "Parallel random numbers: as easy as 1, 2,
 * 3", Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis, 16). Each
 * call maps a 256-bit counter and a 128-bit key to 256 random bits, with no other state. The round constants are those of the
 * Random123 library, and the output matches both that library and numpy.random.Philox.
 *
 * Unsigned 64-bit arithmetic (which wraps modulo 2^64) is done here, rather than in Fortran, which has no unsigned integers and
 * for which signed integer overflow is not defined.
 */

#include <stdint.h>

/* Multipliers and Weyl sequence key increments for Philox4x64. */
#define PHILOX_M4x64_0 UINT64_C(0xD2E7470EE14C6C93)
#define PHILOX_M4x64_1 UINT64_C(0xCA5A826395121157)
#define PHILOX_W64_0   UINT64_C(0x9E3779B97F4A7C15)
#define PHILOX_W64_1   UINT64_C(0xBB67AE8584CAA73B)

static inline void philoxMultiply(uint64_t a, uint64_t b, uint64_t *high, uint64_t *low) {
  /* Return the high and low 64-bit words of the 128-bit product a*b. */
  __uint128_t product = (__uint128_t) a * (__uint128_t) b;
  *high = (uint64_t) (product >> 64);
  *low  = (uint64_t)  product;
}

void philox4x64Block(const uint64_t counter[4], const uint64_t key[2], uint64_t output[4]) {
  /* Apply the ten rounds of Philox4x64 to the given counter and key. */
  uint64_t c0 = counter[0], c1 = counter[1], c2 = counter[2], c3 = counter[3];
  uint64_t k0 = key    [0], k1 = key    [1];
  uint64_t high0, low0, high1, low1;
  int      round;
  for (round = 0; round < 10; round++) {
    if (round > 0) {
      k0 += PHILOX_W64_0;
      k1 += PHILOX_W64_1;
    }
    philoxMultiply(PHILOX_M4x64_0, c0, &high0, &low0);
    philoxMultiply(PHILOX_M4x64_1, c2, &high1, &low1);
    c0 = high1 ^ c1 ^ k0;
    c1 = low1;
    c2 = high0 ^ c3 ^ k1;
    c3 = low0;
  }
  output[0] = c0;
  output[1] = c1;
  output[2] = c2;
  output[3] = c3;
}

void philox4x64Next(uint64_t counter[4], const uint64_t key[2], uint64_t output[4]) {
  /* Increment the 256-bit counter (with carry between words), then generate the block for the new counter value. Incrementing
     before generating matches the sequence produced by numpy.random.Philox. */
  if (++counter[0] == 0)
    if (++counter[1] == 0)
      if (++counter[2] == 0)
        ++counter[3];
  philox4x64Block(counter, key, output);
}
