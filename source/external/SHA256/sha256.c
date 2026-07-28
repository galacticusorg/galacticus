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

/*
  An implementation of the SHA-256 secure hash algorithm, as specified by FIPS 180-4:

    https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.180-4.pdf

  This is used to verify the integrity of content downloaded at run time against a hash pinned in the Galacticus source, for
  cases where the authenticity of the server providing that content can not otherwise be established.
*/

#include "sha256.h"

// Size of the buffer used when reading a file, in bytes.
#define SHA256_READ_BUFFER 65536

// Rotate the 32-bit value `x` right by `n` bits.
#define ROTR(x,n) (((x) >> (n)) | ((x) << (32-(n))))

// The logical functions of FIPS 180-4 section 4.1.2.
#define CH(x,y,z)  (((x) & (y)) ^ (~(x) & (z))            )
#define MAJ(x,y,z) (((x) & (y)) ^ ( (x) & (z)) ^ ((y) & (z)))
#define EP0(x)     (ROTR(x, 2) ^ ROTR(x,13) ^ ROTR(x,22)  )
#define EP1(x)     (ROTR(x, 6) ^ ROTR(x,11) ^ ROTR(x,25)  )
#define SIG0(x)    (ROTR(x, 7) ^ ROTR(x,18) ^ ((x) >>  3) )
#define SIG1(x)    (ROTR(x,17) ^ ROTR(x,19) ^ ((x) >> 10) )

// The first thirty-two bits of the fractional parts of the cube roots of the first sixty-four prime numbers (FIPS 180-4 section
// 4.2.2).
static const uint32_t sha256K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

// Process a single 64-byte block, updating the hash state.
static void sha256Step(uint32_t *state, const uint8_t *block) {
  uint32_t w[64], a, b, c, d, e, f, g, h, temp1, temp2;
  unsigned int i;

  // Prepare the message schedule. The input is interpreted as big-endian, explicitly, so that the result does not depend on the
  // endianness of the host.
  for(i = 0; i < 16; ++i)
    w[i] = ((uint32_t)block[4*i] << 24) | ((uint32_t)block[4*i+1] << 16) | ((uint32_t)block[4*i+2] << 8) | ((uint32_t)block[4*i+3]);
  for(i = 16; i < 64; ++i)
    w[i] = SIG1(w[i-2]) + w[i-7] + SIG0(w[i-15]) + w[i-16];

  a = state[0]; b = state[1]; c = state[2]; d = state[3];
  e = state[4]; f = state[5]; g = state[6]; h = state[7];

  for(i = 0; i < 64; ++i) {
    temp1 = h + EP1(e) + CH(e,f,g) + sha256K[i] + w[i];
    temp2 = EP0(a) + MAJ(a,b,c);
    h = g; g = f; f = e; e = d + temp1;
    d = c; c = b; b = a; a = temp1 + temp2;
  }

  state[0] += a; state[1] += b; state[2] += c; state[3] += d;
  state[4] += e; state[5] += f; state[6] += g; state[7] += h;
}

void sha256Init(SHA256Context *ctx) {
  // The first thirty-two bits of the fractional parts of the square roots of the first eight prime numbers (FIPS 180-4 section
  // 5.3.3).
  ctx->state[0] = 0x6a09e667; ctx->state[1] = 0xbb67ae85;
  ctx->state[2] = 0x3c6ef372; ctx->state[3] = 0xa54ff53a;
  ctx->state[4] = 0x510e527f; ctx->state[5] = 0x9b05688c;
  ctx->state[6] = 0x1f83d9ab; ctx->state[7] = 0x5be0cd19;
  ctx->size        = 0;
  ctx->inputLength = 0;
}

void sha256Update(SHA256Context *ctx, const uint8_t *input, size_t length) {
  size_t i;

  for(i = 0; i < length; ++i) {
    ctx->input[ctx->inputLength] = input[i];
    ++ctx->inputLength;
    ++ctx->size;
    if(ctx->inputLength == 64) {
      sha256Step(ctx->state, ctx->input);
      ctx->inputLength = 0;
    }
  }
}

void sha256Finalize(SHA256Context *ctx, uint8_t *digest) {
  uint64_t     bitCount = ctx->size * 8;
  unsigned int i;

  // Append the "1" bit, then pad with zeros until the length is congruent to 56 modulo 64.
  ctx->input[ctx->inputLength++] = 0x80;
  if(ctx->inputLength > 56) {
    while(ctx->inputLength < 64) ctx->input[ctx->inputLength++] = 0x00;
    sha256Step(ctx->state, ctx->input);
    ctx->inputLength = 0;
  }
  while(ctx->inputLength < 56) ctx->input[ctx->inputLength++] = 0x00;

  // Append the length of the message, in bits, as a big-endian 64-bit integer.
  for(i = 0; i < 8; ++i) ctx->input[56+i] = (uint8_t)(bitCount >> (56 - 8*i));
  sha256Step(ctx->state, ctx->input);

  // Produce the digest, big-endian.
  for(i = 0; i < 8; ++i) {
    digest[4*i  ] = (uint8_t)(ctx->state[i] >> 24);
    digest[4*i+1] = (uint8_t)(ctx->state[i] >> 16);
    digest[4*i+2] = (uint8_t)(ctx->state[i] >>  8);
    digest[4*i+3] = (uint8_t)(ctx->state[i]      );
  }
}

int sha256File(const char *fileName, char *hash) {
  FILE          *file;
  SHA256Context  ctx;
  uint8_t        digest[32];
  uint8_t       *buffer;
  size_t         readCount;
  unsigned int   i;

  file = fopen(fileName, "rb");
  if(file == NULL) return 1;
  buffer = (uint8_t *)malloc(SHA256_READ_BUFFER);
  if(buffer == NULL) {
    fclose(file);
    return 1;
  }
  sha256Init(&ctx);
  while((readCount = fread(buffer, 1, SHA256_READ_BUFFER, file)) > 0)
    sha256Update(&ctx, buffer, readCount);
  if(ferror(file)) {
    free (buffer);
    fclose(file  );
    return 1;
  }
  free (buffer);
  fclose(file  );
  sha256Finalize(&ctx, digest);
  for(i = 0; i < 32; ++i) sprintf(hash+2*i, "%02x", digest[i]);
  hash[64] = '\0';
  return 0;
}
