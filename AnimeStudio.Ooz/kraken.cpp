/*
=== Kraken Decompressor for Windows ===
Copyright (C) 2016, Powzix

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "stdafx.h"

#include <algorithm>
#include <sys/stat.h>

#if defined _WIN32 || defined __CYGWIN__
#ifdef OOZ_DYNAMIC
#ifdef OOZ_BUILD_DLL
#ifdef __GNUC__
#define OOZ_DLL_PUBLIC __attribute__ ((dllexport))
#else
#define OOZ_DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
#endif
#else
#ifdef __GNUC__
#define OOZ_DLL_PUBLIC __attribute__ ((dllimport))
#else
#define OOZ_DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
#endif
#endif
#define OOZ_DLL_LOCAL
#else
#define OOZ_DLL_PUBLIC
#define OOZ_DLL_LOCAL
#endif
#else
#if __GNUC__ >= 4
#define OOZ_DLL_PUBLIC __attribute__ ((visibility ("default")))
#define OOZ_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define OOZ_DLL_PUBLIC
#define OOZ_DLL_LOCAL
#endif
#endif

// c++ is weird
typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned char u8;

// Header in front of each 256k block
typedef struct KrakenHeader {
  int decoder_type; // Type of decoder used, 6 means kraken
  bool restart_decoder; // Whether to restart the decoder
  bool uncompressed; // Whether this block is uncompressed
  bool use_checksums; // Whether this block uses checksums.
  // bool unused_;
} KrakenHeader;

// Kraken decompression happens in two phases, first one decodes
// all the literals and copy lengths using huffman and second
// phase runs the copy loop. This holds the tables needed by stage 2.
typedef struct KrakenLzTable {
  // Stream of (literal, match) pairs. The flag byte contains
  // the length of the match, the length of the literal and whether
  // to use a recent offset.
  byte *cmd_stream;
  int cmd_stream_size;

  // Holds the actual distances in case we're not using a recent offset.
  int *offs_stream;
  int offs_stream_size;

  // Holds the sequence of literals. All literal copying happens from here.
  byte *lit_stream;
  int lit_stream_size;

  // Holds the lengths that do not fit in the flag stream. Both literal
  // lengths and match length are stored in the same array.
  int *len_stream;
  int len_stream_size;
} KrakenLzTable;


// Mermaid/Selkie decompression also happens in two phases, just like in Kraken,
// but the match copier works differently.
// Both Mermaid and Selkie use the same on-disk format, only the compressor
// differs.
typedef struct MermaidLzTable {
  // Flag stream. Format of flags:
  // Read flagbyte from |cmd_stream|
  // If flagbyte >= 24:
  //   flagbyte & 0x80 == 0 : Read from |off16_stream| into |recent_offs|.
  //                   != 0 : Don't read offset.
  //   flagbyte & 7 = Number of literals to copy first from |lit_stream|.
  //   (flagbyte >> 3) & 0xF = Number of bytes to copy from |recent_offs|.
  //
  //  If flagbyte == 0 :
  //    Read byte L from |length_stream|
  //    If L > 251: L += 4 * Read word from |length_stream|
  //    L += 64
  //    Copy L bytes from |lit_stream|.
  //
  //  If flagbyte == 1 :
  //    Read byte L from |length_stream|
  //    If L > 251: L += 4 * Read word from |length_stream|
  //    L += 91
  //    Copy L bytes from match pointed by next offset from |off16_stream|
  //
  //  If flagbyte == 2 :
  //    Read byte L from |length_stream|
  //    If L > 251: L += 4 * Read word from |length_stream|
  //    L += 29
  //    Copy L bytes from match pointed by next offset from |off32_stream|, 
  //    relative to start of block.
  //    Then prefetch |off32_stream[3]|
  //
  //  If flagbyte > 2:
  //    L = flagbyte + 5
  //    Copy L bytes from match pointed by next offset from |off32_stream|,
  //    relative to start of block.
  //    Then prefetch |off32_stream[3]|
  byte *cmd_stream, *cmd_stream_end;
  
  const byte *length_stream; // Length stream
  byte *lit_stream, *lit_stream_end; // Literal stream
  uint16 *off16_stream, *off16_stream_end; // Near offsets
  uint32 *off32_stream, *off32_stream_end; // Far offsets for current chunk
  uint32 *off32_stream_1, *off32_stream_2; // Holds the offsets for the two chunks
  uint32 off32_size_1, off32_size_2;
  uint32 cmd_stream_2_offs, cmd_stream_2_offs_end; // Flag offsets for next 64k chunk.
} MermaidLzTable;

namespace {
typedef struct KrakenDecoder {
  u32 src_len, dst_len, offset, scratch_size;
  const byte *src;
  byte *dst;
  byte *dst1;
  byte *scratch; // Pointer to a 256k buffer that holds the intermediate state in between decode phase 1 and 2.

  KrakenHeader hdr;
} KrakenDecoder;

bool step_ret(KrakenDecoder *dec, const u64 src_used, const u64 dst_used) {
  dec->src += src_used;
  dec->src_len -= src_used;
  dec->dst1 += dst_used;
  dec->offset += dst_used;
  dec->dst_len -= dst_used;
  return src_used > 0;
}
}

typedef struct BitReader {
  const byte *p, *p_end; // |p| holds the current byte and |p_end| the end of the buffer.
  uint32 bits; // Bits accumulated so far
  int bitpos; // Next byte will end up in the |bitpos| position in |bits|.
} BitReader;

struct HuffRevLut {
  uint8 bits2len[2048];
  uint8 bits2sym[2048];
};

typedef struct HuffReader {
  // Array to hold the output of the huffman read array operation
  byte *output, *output_end;
  // We decode three parallel streams, two forwards, |src| and |src_mid|
  // while |src_end| is decoded backwards. 
  const byte *src, *src_mid, *src_end, *src_mid_org;
  int src_bitpos, src_mid_bitpos, src_end_bitpos;
  uint32 src_bits, src_mid_bits, src_end_bits;
} HuffReader;

inline size_t Max(const size_t a, const size_t b) { return a > b ? a : b; }
inline size_t Min(const size_t a, const size_t b) { return a < b ? a : b; }

#define ALIGN_POINTER(p, align) ((uint8*)(((uintptr_t)(p) + (align - 1)) & ~(align - 1)))

struct HuffRange;

int Kraken_DecodeBytes(byte **output, const byte *src, const byte *src_end, int *decoded_size, size_t output_size, bool force_memmove, uint8 *scratch, const uint8 *scratch_end);
int Kraken_GetBlockSize(const uint8 *src, const uint8 *src_end, int *dest_size, int dest_capacity);
int Huff_ConvertToRanges(HuffRange *range, int num_symbols, int P, const uint8 *symlen, BitReader *bits);

// Allocate memory with a specific alignment
void *MallocAligned(const size_t size, const size_t alignment = 0xf) {
  void *x = malloc(size + alignment + sizeof(void*)), *x_org = x;
  if (x) {
    x = (void*)((intptr_t)x + alignment + sizeof(void*) & ~alignment);
    static_cast<void **>(x)[-1] = x_org;
  }
  return x;
}

// Free memory allocated through |MallocAligned|
void FreeAligned(void *p) {
  free(static_cast<void **>(p)[-1]);
}

uint32 BSR(const uint32 x) {
  unsigned long index;
  _BitScanReverse(&index, x);
  return index;
}

uint32 BSF(const uint32 x) {
  unsigned long index;
  _BitScanForward(&index, x);
  return index;
}
int CountLeadingZeros(const uint32 bits) {
  unsigned long x;
  _BitScanReverse(&x, bits);
  return 31 - x;
}

// Read more bytes to make sure we always have at least 24 bits in |bits|.
void BitReader_Refill(BitReader *bits) {
  assert(bits->bitpos <= 24);
  while (bits->bitpos > 0) {
    bits->bits |= (bits->p < bits->p_end ? *bits->p : 0) << bits->bitpos;
    bits->bitpos -= 8;
    bits->p++;
  }
}

// Read more bytes to make sure we always have at least 24 bits in |bits|,
// used when reading backwards.
void BitReader_RefillBackwards(BitReader *bits) {
  assert(bits->bitpos <= 24);
  while (bits->bitpos > 0) {
    bits->p--;
    bits->bits |= (bits->p >= bits->p_end ? *bits->p : 0) << bits->bitpos;
    bits->bitpos -= 8;
  }
}

// Refill bits then read a single bit.
int BitReader_ReadBit(BitReader *bits) {
  BitReader_Refill(bits);
  const int r = bits->bits >> 31;
  bits->bits <<= 1;
  bits->bitpos += 1;
  return r;
}

int BitReader_ReadBitNoRefill(BitReader *bits) {
  const int r = bits->bits >> 31;
  bits->bits <<= 1;
  bits->bitpos += 1;
  return r;
}

// Read |n| bits without refilling.
int BitReader_ReadBitsNoRefill(BitReader *bits, const int n) {
  const int r = bits->bits >> (32 - n);
  bits->bits <<= n;
  bits->bitpos += n;
  return r;
}

// Read |n| bits without refilling, n may be zero.
int BitReader_ReadBitsNoRefillZero(BitReader *bits, const int n) {
  const int r = bits->bits >> 1 >> (31 - n);
  bits->bits <<= n;
  bits->bitpos += n;
  return r;
}

uint32 BitReader_ReadMoreThan24Bits(BitReader *bits, const int n) {
  uint32 rv;
  if (n <= 24) {
    rv = BitReader_ReadBitsNoRefillZero(bits, n);
  } else {
    rv = BitReader_ReadBitsNoRefill(bits, 24) << (n - 24);
    BitReader_Refill(bits);
    rv += BitReader_ReadBitsNoRefill(bits, n - 24);
  }
  BitReader_Refill(bits);
  return rv;
}

uint32 BitReader_ReadMoreThan24BitsB(BitReader *bits, const int n) {
  uint32 rv;
  if (n <= 24) {
    rv = BitReader_ReadBitsNoRefillZero(bits, n);
  } else {
    rv = BitReader_ReadBitsNoRefill(bits, 24) << (n - 24);
    BitReader_RefillBackwards(bits);
    rv += BitReader_ReadBitsNoRefill(bits, n - 24);
  }
  BitReader_RefillBackwards(bits);
  return rv;
}

// Reads a gamma value.
// Assumes bitreader is already filled with at least 23 bits
int BitReader_ReadGamma(BitReader *bits) {
  const int n = 2 * (bits->bits != 0 ? CountLeadingZeros(bits->bits) : 32) + 2;
  assert(n < 24);
  bits->bitpos += n;
  const int r = bits->bits >> (32 - n);
  bits->bits <<= n;
  return r - 2;
}

// Reads a gamma value with |forced| number of forced bits.
int BitReader_ReadGammaX(BitReader *bits, const int forced) {
  if (bits->bits != 0) {
    const int lz = CountLeadingZeros(bits->bits);
    assert(lz < 24);
    const int r = (bits->bits >> (31 - lz - forced)) + ((lz - 1) << forced);
    bits->bits <<= lz + forced + 1;
    bits->bitpos += lz + forced + 1;
    return r;
  }
  return 0;
}

// Reads an offset code parametrized by |v|.
uint32 BitReader_ReadDistance(BitReader *bits, const uint32 v) {
  uint32 w, m, n, rv;
  if (v < 0xF0) {
    n = (v >> 4) + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = ((w & m) << 4) + (v & 0xF) - 248;
  } else {
    n = v - 0xF0 + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = 8322816 + ((w & m) << 12);
    BitReader_Refill(bits);
    rv += bits->bits >> 20;
    bits->bitpos += 12;
    bits->bits <<= 12;
  }
  BitReader_Refill(bits);
  return rv;
}


// Reads a offset code parametrized by |v|, backwards.
uint32 BitReader_ReadDistanceB(BitReader *bits, const uint32 v) {
  uint32 w, m, n, rv;
  if (v < 0xF0) {
    n = (v >> 4) + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = ((w & m) << 4) + (v & 0xF) - 248;
  } else {
    n = v - 0xF0 + 4;
    w = _rotl(bits->bits | 1, n);
    bits->bitpos += n;
    m = (2 << n) - 1;
    bits->bits = w & ~m;
    rv = 8322816 + ((w & m) << 12);
    BitReader_RefillBackwards(bits);
    rv += bits->bits >> (32 - 12);
    bits->bitpos += 12;
    bits->bits <<= 12;
  }
  BitReader_RefillBackwards(bits);
  return rv;
}

// Reads a length code.
bool BitReader_ReadLength(BitReader *bits, uint32 *v) {
  unsigned long bitresult;
  int n;
  uint32 rv;
  _BitScanReverse(&bitresult, bits->bits);
  n = 31 - bitresult;
  if (n > 12) return false;
  bits->bitpos += n;
  bits->bits <<= n;
  BitReader_Refill(bits);
  n += 7;
  bits->bitpos += n;
  rv = (bits->bits >> (32 - n)) - 64;
  bits->bits <<= n;
  *v = rv;
  BitReader_Refill(bits);
  return true;
}

// Reads a length code, backwards.
bool BitReader_ReadLengthB(BitReader *bits, uint32 *v) {
  unsigned long bitresult;
  int n;
  uint32 rv;
  _BitScanReverse(&bitresult, bits->bits);
  n = 31 - bitresult;
  if (n > 12) return false;
  bits->bitpos += n;
  bits->bits <<= n;
  BitReader_RefillBackwards(bits);
  n += 7;
  bits->bitpos += n;
  rv = (bits->bits >> (32 - n)) - 64;
  bits->bits <<= n;
  *v = rv;
  BitReader_RefillBackwards(bits);
  return true;
}

int Log2RoundUp(const uint32 v) {
  if (v > 1) {
    unsigned long idx;
    _BitScanReverse(&idx, v - 1);
    return idx + 1;
  } else {
    return 0;
  }
}

// #define ALIGN_16(x) (((x)+15)&~15)
#define COPY_64(d, s)             \
    {                             \
        uint64_t tmpVal;          \
        memcpy(&tmpVal, (s), 8);  \
        memcpy((d), &tmpVal, 8);  \
    }
#define COPY_64_BYTES(d, s) {                                                 \
        simde_mm_storeu_si128((simde__m128i*)d + 0, simde_mm_loadu_si128((simde__m128i*)s + 0));  \
        simde_mm_storeu_si128((simde__m128i*)d + 1, simde_mm_loadu_si128((simde__m128i*)s + 1));  \
        simde_mm_storeu_si128((simde__m128i*)d + 2, simde_mm_loadu_si128((simde__m128i*)s + 2));  \
        simde_mm_storeu_si128((simde__m128i*)d + 3, simde_mm_loadu_si128((simde__m128i*)s + 3));  \
}
#define COPY_64_ADD(d, s, t) simde_mm_storel_epi64((simde__m128i *)(d), simde_mm_add_epi8(simde_mm_loadl_epi64((simde__m128i *)(s)), simde_mm_loadl_epi64((simde__m128i *)(t))))

static KrakenDecoder *Kraken_Create(const byte *src, const u32 src_len, byte *dst, const u32 dst_len) {
  constexpr size_t scratch_size = 0x6c000;
  const auto dec = static_cast<KrakenDecoder *>(MallocAligned(sizeof(KrakenDecoder) + scratch_size, 16 - 1));
  memset(dec, 0, sizeof(KrakenDecoder));
  dec->src = src;
  dec->src_len = src_len;
  dec->dst = dst;
  dec->dst1 = dst;
  dec->offset = 0;
  dec->dst_len = dst_len;
  dec->scratch_size = scratch_size;
  dec->scratch = (byte*)(dec + 1);
  return dec;
}

void Kraken_Destroy(KrakenDecoder *kraken) {
  FreeAligned(kraken);
}

const byte *LZNA_ParseWholeMatchInfo(const byte *p, uint32 *dist) {
  const uint32 v = _byteswap_ushort(*(uint16*)p);

  if (v < 0x8000) {
    uint32 x = 0, b, pos = 0;
    for (;;) {
      b = p[2];
      p += 1;
      if (b & 0x80)
        break;
      x += (b + 0x80) << pos;
      pos += 7;
    }
    x += (b - 0x80) << pos;
    *dist = 0x8000 + v + (x << 15) + 1;
    return p + 2;
  }
  *dist = v - 0x8000 + 1;
  return p + 2;
}

// Rearranges elements in the input array so that bits in the index
// get flipped.
static void ReverseBitsArray2048(const byte *input, byte *output) {
  static const uint8 offsets[32] = {
    0,    0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8
  };
  simde__m128i t0, t1, t2, t3, s0, s1, s2, s3;
  int i, j;
  for(i = 0; i != 32; i++) {
    j = offsets[i];
    t0 = simde_mm_unpacklo_epi8(simde_mm_loadl_epi64((const simde__m128i *)&input[j]), 
                                simde_mm_loadl_epi64((const simde__m128i *)&input[j + 256]));
    t1 = simde_mm_unpacklo_epi8(simde_mm_loadl_epi64((const simde__m128i *)&input[j + 512]),
                                simde_mm_loadl_epi64((const simde__m128i *)&input[j + 768]));
    t2 = simde_mm_unpacklo_epi8(simde_mm_loadl_epi64((const simde__m128i *)&input[j + 1024]),
                                simde_mm_loadl_epi64((const simde__m128i *)&input[j + 1280]));
    t3 = simde_mm_unpacklo_epi8(simde_mm_loadl_epi64((const simde__m128i *)&input[j + 1536]),
                                simde_mm_loadl_epi64((const simde__m128i *)&input[j + 1792]));

    s0 = simde_mm_unpacklo_epi8(t0, t1);
    s1 = simde_mm_unpacklo_epi8(t2, t3);
    s2 = simde_mm_unpackhi_epi8(t0, t1);
    s3 = simde_mm_unpackhi_epi8(t2, t3);

    t0 = simde_mm_unpacklo_epi8(s0, s1);
    t1 = simde_mm_unpacklo_epi8(s2, s3);
    t2 = simde_mm_unpackhi_epi8(s0, s1);
    t3 = simde_mm_unpackhi_epi8(s2, s3);

    simde_mm_storel_epi64((simde__m128i *)&output[0], t0);
    simde_mm_storeh_pi((simde__m64*)&output[1024], simde_mm_castsi128_ps(t0));
    simde_mm_storel_epi64((simde__m128i *)&output[256], t1);
    simde_mm_storeh_pi((simde__m64*)&output[1280], simde_mm_castsi128_ps(t1));
    simde_mm_storel_epi64((simde__m128i *)&output[512], t2);
    simde_mm_storeh_pi((simde__m64*)&output[1536], simde_mm_castsi128_ps(t2));
    simde_mm_storel_epi64((simde__m128i *)&output[768], t3);
    simde_mm_storeh_pi((simde__m64*)&output[1792], simde_mm_castsi128_ps(t3));
    output += 8;
  }
}

bool Kraken_DecodeBytesCore(const HuffReader *hr, const HuffRevLut *lut) {
  const byte *src = hr->src;
  uint32 src_bits = hr->src_bits;
  int src_bitpos = hr->src_bitpos;

  const byte *src_mid = hr->src_mid;
  uint32 src_mid_bits = hr->src_mid_bits;
  int src_mid_bitpos = hr->src_mid_bitpos;

  const byte *src_end = hr->src_end;
  uint32 src_end_bits = hr->src_end_bits;
  int src_end_bitpos = hr->src_end_bitpos;

  int k, n;

  byte *dst = hr->output;
  const byte *dst_end = hr->output_end;

  if (src > src_mid)
    return false;

  if (hr->src_end - src_mid >= 4 && dst_end - dst >= 6) {
    dst_end -= 5;
    src_end -= 4;

    while (dst < dst_end && src <= src_mid && src_mid <= src_end) {
      src_bits |= *(uint32*)src << src_bitpos;
      src += (31 - src_bitpos) >> 3;

      src_end_bits |= _byteswap_ulong(*(uint32*)src_end) << src_end_bitpos;
      src_end -= (31 - src_end_bitpos) >> 3;

      src_mid_bits |= *(uint32*)src_mid << src_mid_bitpos;
      src_mid += (31 - src_mid_bitpos) >> 3;

      src_bitpos |= 0x18;
      src_end_bitpos |= 0x18;
      src_mid_bitpos |= 0x18;

      k = src_bits & 0x7FF;
      n = lut->bits2len[k];
      src_bits >>= n;
      src_bitpos -= n;
      dst[0] = lut->bits2sym[k];

      k = src_end_bits & 0x7FF;
      n = lut->bits2len[k];
      src_end_bits >>= n;
      src_end_bitpos -= n;
      dst[1] = lut->bits2sym[k];

      k = src_mid_bits & 0x7FF;
      n = lut->bits2len[k];
      src_mid_bits >>= n;
      src_mid_bitpos -= n;
      dst[2] = lut->bits2sym[k];

      k = src_bits & 0x7FF;
      n = lut->bits2len[k];
      src_bits >>= n;
      src_bitpos -= n;
      dst[3] = lut->bits2sym[k];

      k = src_end_bits & 0x7FF;
      n = lut->bits2len[k];
      src_end_bits >>= n;
      src_end_bitpos -= n;
      dst[4] = lut->bits2sym[k];

      k = src_mid_bits & 0x7FF;
      n = lut->bits2len[k];
      src_mid_bits >>= n;
      src_mid_bitpos -= n;
      dst[5] = lut->bits2sym[k];
      dst += 6;
    }
    dst_end += 5;

    src -= src_bitpos >> 3;
    src_bitpos &= 7;

    src_end += 4 + (src_end_bitpos >> 3);
    src_end_bitpos &= 7;

    src_mid -= src_mid_bitpos >> 3;
    src_mid_bitpos &= 7;
  }
  for(;;) {
    if (dst >= dst_end)
      break;

    if (src_mid - src <= 1) {
      if (src_mid - src == 1)
        src_bits |= *src << src_bitpos;
    } else {
      src_bits |= *(uint16 *)src << src_bitpos;
    }
    k = src_bits & 0x7FF;
    n = lut->bits2len[k];
    src_bitpos -= n;
    src_bits >>= n;
    *dst++ = lut->bits2sym[k];
    src += (7 - src_bitpos) >> 3;
    src_bitpos &= 7;

    if (dst < dst_end) {
      if (src_end - src_mid <= 1) {
        if (src_end - src_mid == 1) {
          src_end_bits |= *src_mid << src_end_bitpos;
          src_mid_bits |= *src_mid << src_mid_bitpos;
        }
      } else {
        const unsigned int v = *(uint16*)(src_end - 2);
        src_end_bits |= (((v >> 8) | (v << 8)) & 0xffff) << src_end_bitpos;
        src_mid_bits |= *(uint16*)src_mid << src_mid_bitpos;
      }
      n = lut->bits2len[src_end_bits & 0x7FF];
      *dst++ = lut->bits2sym[src_end_bits & 0x7FF];
      src_end_bitpos -= n;
      src_end_bits >>= n;
      src_end -= (7 - src_end_bitpos) >> 3;
      src_end_bitpos &= 7;
      if (dst < dst_end) {
        n = lut->bits2len[src_mid_bits & 0x7FF];
        *dst++ = lut->bits2sym[src_mid_bits & 0x7FF];
        src_mid_bitpos -= n;
        src_mid_bits >>= n;
        src_mid += (7 - src_mid_bitpos) >> 3;
        src_mid_bitpos &= 7;
      }
    }
    if (src > src_mid || src_mid > src_end)
      return false;
  }
  if (src != hr->src_mid_org || src_end != src_mid)
    return false;
  return true;
}

int Huff_ReadCodeLengthsOld(BitReader *bits, uint8 *syms, uint32 *code_prefix) {
  if (BitReader_ReadBitNoRefill(bits)) {
    int n, sym = 0, codelen, num_symbols = 0;
    int avg_bits_x4 = 32;
    const int forced_bits = BitReader_ReadBitsNoRefill(bits, 2);

    const uint32 thres_for_valid_gamma_bits = 1 << (31 - (20u >> forced_bits));
    if (BitReader_ReadBit(bits))
      goto SKIP_INITIAL_ZEROS;
    do {
      // Run of zeros
      if (!(bits->bits & 0xff000000))
        return -1;
      sym += BitReader_ReadBitsNoRefill(bits, 2 * (CountLeadingZeros(bits->bits) + 1)) - 2 + 1;
      if (sym >= 256)
        break;
    SKIP_INITIAL_ZEROS:
      BitReader_Refill(bits);
      // Read out the gamma value for the # of symbols
      if (!(bits->bits & 0xff000000))
        return -1;
      n = BitReader_ReadBitsNoRefill(bits, 2 * (CountLeadingZeros(bits->bits) + 1)) - 2 + 1;
      // Overflow?
      if (sym + n > 256)
        return -1;
      BitReader_Refill(bits);
      num_symbols += n;
      do {
        if (bits->bits < thres_for_valid_gamma_bits)
          return -1; // too big gamma value?

        const int lz = CountLeadingZeros(bits->bits);
        const int v = BitReader_ReadBitsNoRefill(bits, lz + forced_bits + 1) + ((lz - 1) << forced_bits);
        codelen = (-(int)(v & 1) ^ (v >> 1)) + ((avg_bits_x4 + 2) >> 2);
        if (codelen < 1 || codelen > 11)
          return -1;
        avg_bits_x4 = codelen + ((3 * avg_bits_x4 + 2) >> 2);
        BitReader_Refill(bits);
        syms[code_prefix[codelen]++] = sym++;
      } while (--n);
    } while (sym != 256);
    return sym == 256 && num_symbols >= 2 ? num_symbols : -1;
  } else {
    // Sparse symbol encoding
    const int num_symbols = BitReader_ReadBitsNoRefill(bits, 8);
    if (num_symbols == 0)
      return -1;
    if (num_symbols == 1) {
      syms[0] = BitReader_ReadBitsNoRefill(bits, 8);
    } else {
      const int codelen_bits = BitReader_ReadBitsNoRefill(bits, 3);
      if (codelen_bits > 4)
        return -1;
      for (int i = 0; i < num_symbols; i++) {
        BitReader_Refill(bits);
        const int sym = BitReader_ReadBitsNoRefill(bits, 8);
        const int codelen = BitReader_ReadBitsNoRefillZero(bits, codelen_bits) + 1;
        if (codelen > 11)
          return -1;
        syms[code_prefix[codelen]++] = sym;
      }
    }
    return num_symbols;
  }
}

int BitReader_ReadFluff(BitReader *bits, const int num_symbols) {
  unsigned long y;

  if (num_symbols == 256)
    return 0;

  int x = 257 - num_symbols;
  if (x > num_symbols)
    x = num_symbols;

  x *= 2;

  _BitScanReverse(&y, x - 1);
  y += 1;

  const uint32 v = bits->bits >> (32 - y);
  const uint32 z = (1 << y) - x;

  if (v >> 1 >= z) {
    bits->bits <<= y;
    bits->bitpos += y;
    return v - z;
  } else {
    bits->bits <<= y - 1;
    bits->bitpos += y - 1;
    return v >> 1;
  }
}

struct BitReader2 {
  const uint8 *p, *p_end;
  uint32 bitpos;
};

static const uint32 kRiceCodeBits2Value[256] = {
  0x80000000, 0x00000007, 0x10000006, 0x00000006, 0x20000005, 0x00000105, 0x10000005, 0x00000005,
  0x30000004, 0x00000204, 0x10000104, 0x00000104, 0x20000004, 0x00010004, 0x10000004, 0x00000004,
  0x40000003, 0x00000303, 0x10000203, 0x00000203, 0x20000103, 0x00010103, 0x10000103, 0x00000103,
  0x30000003, 0x00020003, 0x10010003, 0x00010003, 0x20000003, 0x01000003, 0x10000003, 0x00000003,
  0x50000002, 0x00000402, 0x10000302, 0x00000302, 0x20000202, 0x00010202, 0x10000202, 0x00000202,
  0x30000102, 0x00020102, 0x10010102, 0x00010102, 0x20000102, 0x01000102, 0x10000102, 0x00000102,
  0x40000002, 0x00030002, 0x10020002, 0x00020002, 0x20010002, 0x01010002, 0x10010002, 0x00010002,
  0x30000002, 0x02000002, 0x11000002, 0x01000002, 0x20000002, 0x00000012, 0x10000002, 0x00000002,
  0x60000001, 0x00000501, 0x10000401, 0x00000401, 0x20000301, 0x00010301, 0x10000301, 0x00000301,
  0x30000201, 0x00020201, 0x10010201, 0x00010201, 0x20000201, 0x01000201, 0x10000201, 0x00000201,
  0x40000101, 0x00030101, 0x10020101, 0x00020101, 0x20010101, 0x01010101, 0x10010101, 0x00010101,
  0x30000101, 0x02000101, 0x11000101, 0x01000101, 0x20000101, 0x00000111, 0x10000101, 0x00000101,
  0x50000001, 0x00040001, 0x10030001, 0x00030001, 0x20020001, 0x01020001, 0x10020001, 0x00020001,
  0x30010001, 0x02010001, 0x11010001, 0x01010001, 0x20010001, 0x00010011, 0x10010001, 0x00010001,
  0x40000001, 0x03000001, 0x12000001, 0x02000001, 0x21000001, 0x01000011, 0x11000001, 0x01000001,
  0x30000001, 0x00000021, 0x10000011, 0x00000011, 0x20000001, 0x00001001, 0x10000001, 0x00000001,
  0x70000000, 0x00000600, 0x10000500, 0x00000500, 0x20000400, 0x00010400, 0x10000400, 0x00000400,
  0x30000300, 0x00020300, 0x10010300, 0x00010300, 0x20000300, 0x01000300, 0x10000300, 0x00000300,
  0x40000200, 0x00030200, 0x10020200, 0x00020200, 0x20010200, 0x01010200, 0x10010200, 0x00010200,
  0x30000200, 0x02000200, 0x11000200, 0x01000200, 0x20000200, 0x00000210, 0x10000200, 0x00000200,
  0x50000100, 0x00040100, 0x10030100, 0x00030100, 0x20020100, 0x01020100, 0x10020100, 0x00020100,
  0x30010100, 0x02010100, 0x11010100, 0x01010100, 0x20010100, 0x00010110, 0x10010100, 0x00010100,
  0x40000100, 0x03000100, 0x12000100, 0x02000100, 0x21000100, 0x01000110, 0x11000100, 0x01000100,
  0x30000100, 0x00000120, 0x10000110, 0x00000110, 0x20000100, 0x00001100, 0x10000100, 0x00000100,
  0x60000000, 0x00050000, 0x10040000, 0x00040000, 0x20030000, 0x01030000, 0x10030000, 0x00030000,
  0x30020000, 0x02020000, 0x11020000, 0x01020000, 0x20020000, 0x00020010, 0x10020000, 0x00020000,
  0x40010000, 0x03010000, 0x12010000, 0x02010000, 0x21010000, 0x01010010, 0x11010000, 0x01010000,
  0x30010000, 0x00010020, 0x10010010, 0x00010010, 0x20010000, 0x00011000, 0x10010000, 0x00010000,
  0x50000000, 0x04000000, 0x13000000, 0x03000000, 0x22000000, 0x02000010, 0x12000000, 0x02000000,
  0x31000000, 0x01000020, 0x11000010, 0x01000010, 0x21000000, 0x01001000, 0x11000000, 0x01000000,
  0x40000000, 0x00000030, 0x10000020, 0x00000020, 0x20000010, 0x00001010, 0x10000010, 0x00000010,
  0x30000000, 0x00002000, 0x10001000, 0x00001000, 0x20000000, 0x00100000, 0x10000000, 0x00000000,
};

static const uint8 kRiceCodeBits2Len[256] = {
  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
};


bool DecodeGolombRiceLengths(uint8 *dst, const size_t size, BitReader2 *br) {
  const uint8 *p = br->p, *p_end = br->p_end;
  const uint8 *dst_end = dst + size;
  if (p >= p_end)
    return false;
   
  int count = -static_cast<int>(br->bitpos);
  uint32 v = *p++ & (255 >> br->bitpos);
  for (;;) {
    if (v == 0) {
      count += 8;
    } else {
      const uint32 x = kRiceCodeBits2Value[v];
      *(uint32*)&dst[0] = count + (x & 0x0f0f0f0f);
      *(uint32*)&dst[4] = (x >> 4) & 0x0f0f0f0f;
      dst += kRiceCodeBits2Len[v];
      if (dst >= dst_end)
        break;
      count = x >> 28;
    }
    if (p >= p_end)
      return false;
    v = *p++;
  }
  // went too far, step back
  if (dst > dst_end) {
    int n = dst - dst_end;
    do v &= v - 1; while (--n);
  }
  // step back if byte not finished
  int bitpos = 0;
  if (!(v & 1)) {
    p--;
    unsigned long q;
    _BitScanForward(&q, v);
    bitpos = 8 - q;
  }
  br->p = p;
  br->bitpos = bitpos;
  return true;
}

bool DecodeGolombRiceBits(uint8 *dst, const uint size, const uint bitcount, BitReader2 *br) {
  if (bitcount == 0)
    return true;
  uint8 *dst_end = dst + size;
  const uint8 *p = br->p;
  const int bitpos = br->bitpos;

  const uint bits_required = bitpos + bitcount * size;
  const uint bytes_required = (bits_required + 7) >> 3;
  if (bytes_required > br->p_end - p)
    return false;

  br->p = p + (bits_required >> 3);
  br->bitpos = bits_required & 7;

  // todo. handle r/w outside of range
  const uint64 bak = *(uint64*)dst_end;

  if (bitcount < 2) {
    assert(bitcount == 1);
    do {
      // Read the next byte
      uint64 bits = static_cast<uint8>(_byteswap_ulong(*(uint32 *)p) >> (24 - bitpos));
      p += 1;
      // Expand each bit into each byte of the uint64.
      bits = (bits | (bits << 28)) & 0xF0000000Full;
      bits = (bits | (bits << 14)) & 0x3000300030003ull;
      bits = (bits | (bits <<  7)) & 0x0101010101010101ull;
      *(uint64*)dst = *(uint64*)dst * 2 + _byteswap_uint64(bits);
      dst += 8;
    } while (dst < dst_end);
  } else if (bitcount == 2) {
    do {
      // Read the next 2 bytes
      uint64 bits = static_cast<uint16>(_byteswap_ulong(*(uint32 *)p) >> (16 - bitpos));
      p += 2;
      // Expand each bit into each byte of the uint64.
      bits = (bits | (bits << 24)) & 0xFF000000FFull;
      bits = (bits | (bits << 12)) & 0xF000F000F000Full;
      bits = (bits | (bits << 6)) & 0x0303030303030303ull;
      *(uint64*)dst = *(uint64*)dst * 4 + _byteswap_uint64(bits);
      dst += 8;
    } while (dst < dst_end);

  } else {
    assert(bitcount == 3);
    do {
      // Read the next 3 bytes
      uint64 bits = (_byteswap_ulong(*(uint32*)p) >> (8 - bitpos)) & 0xffffff;
      p += 3;
      // Expand each bit into each byte of the uint64.
      bits = (bits | (bits << 20)) & 0xFFF00000FFFull;
      bits = (bits | (bits << 10)) & 0x3F003F003F003Full;
      bits = (bits | (bits << 5)) & 0x0707070707070707ull;
      *(uint64*)dst = *(uint64*)dst * 8 + _byteswap_uint64(bits);
      dst += 8;
    } while (dst < dst_end);
  }
  *(uint64*)dst_end = bak;
  return true;
}

struct HuffRange {
  uint16 symbol;
  uint16 num;
};

int Huff_ConvertToRanges(HuffRange *range, const int num_symbols, const int P, const uint8 *symlen, BitReader *bits) {
  int num_ranges = P >> 1, v, sym_idx = 0;

  // Start with space?
  if (P & 1) {
    BitReader_Refill(bits);
    v = *symlen++;
    if (v >= 8)
      return -1;
    sym_idx = BitReader_ReadBitsNoRefill(bits, v + 1) + (1 << (v + 1)) - 1;
  }
  int syms_used = 0;

  for (int i = 0; i < num_ranges; i++) {
    BitReader_Refill(bits);
    v = symlen[0];
    if (v >= 9)
      return -1;
    const int num = BitReader_ReadBitsNoRefillZero(bits, v) + (1 << v);
    v = symlen[1];
    if (v >= 8)
      return -1;
    const int space = BitReader_ReadBitsNoRefill(bits, v + 1) + (1 << (v + 1)) - 1;
    range[i].symbol = sym_idx;
    range[i].num = num;
    syms_used += num;
    sym_idx += num + space;
    symlen += 2;
  }

  if (sym_idx >= 256 || syms_used >= num_symbols || sym_idx + num_symbols - syms_used > 256)
    return -1;

  range[num_ranges].symbol = sym_idx;
  range[num_ranges].num = num_symbols - syms_used;

  return num_ranges + 1;
}

int Huff_ReadCodeLengthsNew(BitReader *bits, uint8 *syms, uint32 *code_prefix) {
  const int forced_bits = BitReader_ReadBitsNoRefill(bits, 2);
  const int num_symbols = BitReader_ReadBitsNoRefill(bits, 8) + 1;
  const int fluff = BitReader_ReadFluff(bits, num_symbols);

  uint8 code_len[512];
  BitReader2 br2;
  br2.bitpos = bits->bitpos - 24 & 7;
  br2.p_end = bits->p_end;
  br2.p = bits->p - static_cast<unsigned>((24 - bits->bitpos + 7) >> 3);

  if (!DecodeGolombRiceLengths(code_len, num_symbols + fluff, &br2)) return -1;
  memset(code_len + (num_symbols + fluff), 0, 16);
  if (!DecodeGolombRiceBits(code_len, num_symbols, forced_bits, &br2)) return -1;
   
  // Reset the bits decoder.
  bits->bitpos = 24;
  bits->p = br2.p;
  bits->bits = 0;
  BitReader_Refill(bits);
  bits->bits <<= br2.bitpos;
  bits->bitpos += br2.bitpos;

  uint running_sum = 0x1e;
  for (int i = 0; i < num_symbols; i++) {
    int v = code_len[i];
    v = -(v & 1) ^ (v >> 1);
    code_len[i] = v + (running_sum >> 2) + 1;
    if (code_len[i] < 1 || code_len[i] > 11) return -1;
    running_sum += v;
  }

  HuffRange range[128];
  const int ranges = Huff_ConvertToRanges(range, num_symbols, fluff, &code_len[num_symbols], bits);
  if (ranges <= 0) return -1;
  uint8 *cp = code_len;
  for (int i = 0; i < ranges; i++) {
    int sym = range[i].symbol;
    int n = range[i].num;
    do {
      syms[code_prefix[*cp++]++] = sym++;
    } while (--n);
  }

  return num_symbols;
}

struct NewHuffLut {
  uint8 bits2len[2048 + 16]; // Mapping that maps a bit pattern to a code length.
  uint8 bits2sym[2048 + 16]; // Mapping that maps a bit pattern to a symbol.
};

// May overflow 16 bytes past the end
void FillByteOverflow16(uint8 *dst, const uint8 v, const size_t n) {
  memset(dst, v, n);
}

bool Huff_MakeLut(const uint32 *prefix_org, const uint32 *prefix_cur, NewHuffLut *hufflut, const uint8 *syms) {
  uint32 currslot = 0;
  for(uint32 i = 1; i < 11; i++) {
    const uint32 start = prefix_org[i];
    const uint32 count = prefix_cur[i] - start;
    if (count) {
      const uint32 stepsize = 1 << (11 - i);
      const uint32 num_to_set = count << (11 - i);
      if (currslot + num_to_set > 2048) return false;
      FillByteOverflow16(&hufflut->bits2len[currslot], i, num_to_set);

      uint8 *p = &hufflut->bits2sym[currslot];
      for (uint32 j = 0; j != count; j++, p += stepsize)
        FillByteOverflow16(p, syms[start + j], stepsize);
      currslot += num_to_set;
    }
  }
  if (prefix_cur[11] - prefix_org[11] != 0) {
    const uint32 num_to_set = prefix_cur[11] - prefix_org[11];
    if (currslot + num_to_set > 2048) return false;
    FillByteOverflow16(&hufflut->bits2len[currslot], 11, num_to_set);
    memcpy(&hufflut->bits2sym[currslot], &syms[prefix_org[11]], num_to_set);
    currslot += num_to_set;
  }
  return currslot == 2048;
}

int Kraken_DecodeBytes_Type12(const byte *src, const size_t src_size, byte *output, const int output_size, const int type) {
  int half_output_size;
  uint32 split_left, split_mid, split_right;
  const byte *src_mid;
  NewHuffLut huff_lut;
  HuffReader hr;
  HuffRevLut rev_lut;
  const uint8 *src_end = src + src_size;

  BitReader bits;
  bits.bitpos = 24;
  bits.bits = 0;
  bits.p = src;
  bits.p_end = src_end;
  BitReader_Refill(&bits);

  static const uint32 code_prefix_org[12] = { 0x0, 0x0, 0x2, 0x6, 0xE, 0x1E, 0x3E, 0x7E, 0xFE, 0x1FE, 0x2FE, 0x3FE };
  uint32 code_prefix[12] = { 0x0, 0x0, 0x2, 0x6, 0xE, 0x1E, 0x3E, 0x7E, 0xFE, 0x1FE, 0x2FE, 0x3FE };
  uint8 syms[1280];
  int num_syms;
  if (!BitReader_ReadBitNoRefill(&bits)) {
    num_syms = Huff_ReadCodeLengthsOld(&bits, syms, code_prefix);
  } else if (!BitReader_ReadBitNoRefill(&bits)) {
    num_syms = Huff_ReadCodeLengthsNew(&bits, syms, code_prefix);
  } else {
    return -1;
  }

  if (num_syms < 1) return -1;
  src = bits.p - (24 - bits.bitpos) / 8;

  if (num_syms == 1) {
    memset(output, syms[0], output_size);
    return src - src_end;
  }
  
  if (!Huff_MakeLut(code_prefix_org, code_prefix, &huff_lut, syms)) return -1;
  ReverseBitsArray2048(huff_lut.bits2len, rev_lut.bits2len);
  ReverseBitsArray2048(huff_lut.bits2sym, rev_lut.bits2sym);

  if (type == 1) {
    if (src + 3 > src_end) return -1;
    split_mid = *(uint16*)src;
    src += 2;
    hr.output = output;
    hr.output_end = output + output_size;
    hr.src = src;
    hr.src_end = src_end;
    hr.src_mid_org = hr.src_mid = src + split_mid;
    hr.src_bitpos = 0;
    hr.src_bits = 0;
    hr.src_mid_bitpos = 0;
    hr.src_mid_bits = 0;
    hr.src_end_bitpos = 0;
    hr.src_end_bits = 0;
    if (!Kraken_DecodeBytesCore(&hr, &rev_lut)) return -1;
  } else {
    if (src + 6 > src_end) return -1;

    half_output_size = (output_size + 1) >> 1;
    split_mid = *(uint32*)src & 0xFFFFFF;
    src += 3;
    if (split_mid > src_end - src) return -1;
    src_mid = src + split_mid;
    split_left = *(uint16*)src;
    src += 2;
    if (src_mid - src < split_left + 2 || src_end - src_mid < 3) return -1;
    split_right = *(uint16*)src_mid;
    if (src_end - (src_mid + 2) < split_right + 2) return -1;

    hr.output = output;
    hr.output_end = output + half_output_size;
    hr.src = src;
    hr.src_end = src_mid;
    hr.src_mid_org = hr.src_mid = src + split_left;
    hr.src_bitpos = 0;
    hr.src_bits = 0;
    hr.src_mid_bitpos = 0;
    hr.src_mid_bits = 0;
    hr.src_end_bitpos = 0;
    hr.src_end_bits = 0;
    if (!Kraken_DecodeBytesCore(&hr, &rev_lut)) return -1;

    hr.output = output + half_output_size;
    hr.output_end = output + output_size;
    hr.src = src_mid + 2;
    hr.src_end = src_end;
    hr.src_mid_org = hr.src_mid = src_mid + 2 + split_right;
    hr.src_bitpos = 0;
    hr.src_bits = 0;
    hr.src_mid_bitpos = 0;
    hr.src_mid_bits = 0;
    hr.src_end_bitpos = 0;
    hr.src_end_bits = 0;
    if (!Kraken_DecodeBytesCore(&hr, &rev_lut)) return -1;
  }
  return static_cast<int>(src_size);
}

static uint32 bitmasks[32] = {
  0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff,
  0x1ff, 0x3ff, 0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff,
  0x1ffff, 0x3ffff, 0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff,
  0xffffff, 0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff
};

int Kraken_DecodeMultiArray(const uint8 *src, const uint8 *src_end,
                            uint8 *dst, uint8 *dst_end,
                            uint8 **array_data, int *array_lens, int array_count,
                            int *total_size_out, bool force_memmove, uint8 *scratch, const uint8 *scratch_end) {
  const uint8 *src_org = src;
  if (src_end - src < 4) return -1;
  int decoded_size;
  int num_arrays_in_file = *src++;
  if (!(num_arrays_in_file & 0x80)) return -1;
  num_arrays_in_file &= 0x3f;

  if (dst == scratch) {
    // todo: ensure scratch space first?
    scratch += (scratch_end - scratch - 0xc000) >> 1;
    dst_end = scratch;
  }

  int total_size = 0;
  if (num_arrays_in_file == 0) {
    for (int i = 0; i < array_count; i++) {
      uint8 *chunk_dst = dst;
      int dec = Kraken_DecodeBytes(&chunk_dst, src, src_end, &decoded_size, dst_end - dst, force_memmove, scratch, scratch_end);
      if (dec < 0) return -1;
      dst += decoded_size;
      array_lens[i] = decoded_size;
      array_data[i] = chunk_dst;
      src += dec;
      total_size += decoded_size;
    }
    *total_size_out = total_size;
    return src - src_org; // not supported yet
  }

  uint8 *entropy_array_data[63];
  uint32 entropy_array_size[63];

  // First loop just decodes everything to scratch
  uint8 *scratch_cur = scratch;

  for(int i = 0; i < num_arrays_in_file; i++) {
    uint8 *chunk_dst = scratch_cur;
    int dec = Kraken_DecodeBytes(&chunk_dst, src, src_end, &decoded_size, scratch_end - scratch_cur, force_memmove, scratch_cur, scratch_end);
    if (dec < 0) return -1;
    entropy_array_data[i] = chunk_dst;
    entropy_array_size[i] = decoded_size;
    scratch_cur += decoded_size;
    total_size += decoded_size;
    src += dec;
  }
  *total_size_out = total_size;

  if (src_end - src < 3) return -1;
  int Q = *(uint16*)src;
  src += 2;
  int out_size;
  if (Kraken_GetBlockSize(src, src_end, &out_size, total_size) < 0) return -1;
  int num_indexes = out_size;

  int num_lens = num_indexes - array_count;
  if (num_lens < 1) return -1;
  if (scratch_end - scratch_cur < num_indexes) return -1;
  uint8 *interval_lenlog2 = scratch_cur;
  scratch_cur += num_indexes;

  if (scratch_end - scratch_cur < num_indexes) return -1;
  uint8 *interval_indexes = scratch_cur;
  scratch_cur += num_indexes;
 
  if (Q & 0x8000) {
    int size_out;
    int n = Kraken_DecodeBytes(&interval_indexes, src, src_end, &size_out, num_indexes, true, scratch_cur, scratch_end);
    if (n < 0 || size_out != num_indexes) return -1;
    src += n;

    for (int i = 0; i < num_indexes; i++) {
      int t = interval_indexes[i];
      interval_lenlog2[i] = t >> 4;
      interval_indexes[i] = t & 0xF;
    }

    num_lens = num_indexes;
  } else {
    int lenlog2_chunksize = num_indexes - array_count;
    int size_out;
    int n = Kraken_DecodeBytes(&interval_indexes, src, src_end, &size_out, num_indexes, false, scratch_cur, scratch_end);
    if (n < 0 || size_out != num_indexes) return -1;
    src += n;

    n = Kraken_DecodeBytes(&interval_lenlog2, src, src_end, &size_out, lenlog2_chunksize, false, scratch_cur, scratch_end);
    if (n < 0 || size_out != lenlog2_chunksize) return -1;
    src += n;

    for (int i = 0; i < lenlog2_chunksize; i++)
      if (interval_lenlog2[i] > 16) return -1;
  }

  if (scratch_end - scratch_cur < 4) return -1;
  scratch_cur = ALIGN_POINTER(scratch_cur, 4);
  if (scratch_end - scratch_cur < num_lens * 4) return -1;
  uint32 *decoded_intervals = (uint32*)scratch_cur;
  int varbits_complen = Q & 0x3FFF;
  if (src_end - src < varbits_complen) return -1;
  const uint8 *f = src;
  uint32 bits_f = 0;
  int bitpos_f = 24;
  const uint8 *src_end_actual = src + varbits_complen;
  const uint8 *b = src_end_actual;
  uint32 bits_b = 0;
  int bitpos_b = 24;

  int i;
  for (i = 0; i + 2 <= num_lens; i += 2) {
    bits_f |= _byteswap_ulong(*(uint32*)f) >> (24 - bitpos_f);
    f += (bitpos_f + 7) >> 3;

    bits_b |= ((uint32*)b)[-1] >> (24 - bitpos_b);
    b -= (bitpos_b + 7) >> 3;

    int numbits_f = interval_lenlog2[i + 0];
    int numbits_b = interval_lenlog2[i + 1];

    bits_f = _rotl(bits_f | 1, numbits_f);
    bitpos_f += numbits_f - 8 * ((bitpos_f + 7) >> 3);

    bits_b = _rotl(bits_b | 1, numbits_b);
    bitpos_b += numbits_b - 8 * ((bitpos_b + 7) >> 3);

    int value_f = bits_f & bitmasks[numbits_f];
    bits_f &= ~bitmasks[numbits_f];

    int value_b = bits_b & bitmasks[numbits_b];
    bits_b &= ~bitmasks[numbits_b];

    decoded_intervals[i + 0] = value_f;
    decoded_intervals[i + 1] = value_b;
  }

  // read final one since above loop reads 2
  if (i < num_lens) {
    bits_f |= _byteswap_ulong(*(uint32*)f) >> (24 - bitpos_f);
    int numbits_f = interval_lenlog2[i];
    bits_f = _rotl(bits_f | 1, numbits_f);
    int value_f = bits_f & bitmasks[numbits_f];
    decoded_intervals[i + 0] = value_f;
  }

  if (interval_indexes[num_indexes - 1]) return -1;
  int indi = 0, leni = 0, source;
  int increment_leni = (Q & 0x8000) != 0;

  for(int arri = 0; arri < array_count; arri++) {
    array_data[arri] = dst;
    if (indi >= num_indexes) return -1;

    while ((source = interval_indexes[indi++]) != 0) {
      if (source > num_arrays_in_file) return -1;
      if (leni >= num_lens) return -1;
      int cur_len = decoded_intervals[leni++];
      int bytes_left = entropy_array_size[source - 1];
      if (cur_len > bytes_left || cur_len > dst_end - dst) return -1;
      uint8 *blksrc = entropy_array_data[source - 1];
      entropy_array_size[source - 1] -= cur_len;
      entropy_array_data[source - 1] += cur_len;
      uint8 *dstx = dst;
      dst += cur_len;
      memcpy(dstx, blksrc, cur_len);
    }
    leni += increment_leni;
    array_lens[arri] = dst - array_data[arri];
  }

  if (indi != num_indexes || leni != num_lens) return -1;
  for (int i = 0; i < num_arrays_in_file; i++) {
    if (entropy_array_size[i]) return -1;
  }
  return src_end_actual - src_org;
}

int Krak_DecodeRecursive(const byte *src, const size_t src_size, byte *output, const int output_size, uint8 *scratch, const uint8 *scratch_end) {
  const uint8 *src_org = src;
  byte *output_end = output + output_size;
  const byte *src_end = src + src_size;

  if (src_size < 6)
    return -1;

  int n = src[0] & 0x7f;
  if (n < 2)
    return -1;

  if (!(src[0] & 0x80)) {
    src++;
    do {
      int decoded_size;
      const int dec = Kraken_DecodeBytes(&output, src, src_end, &decoded_size, output_end - output, true, scratch, scratch_end);
      if (dec < 0)
        return -1;
      output += decoded_size;
      src += dec;
    } while (--n);
    if (output != output_end)
      return -1;
    return src - src_org;
  } else {
    uint8 *array_data;
    int array_len, decoded_size;
    const int dec = Kraken_DecodeMultiArray(src, src_end, output, output_end, &array_data, &array_len, 1, &decoded_size, true, scratch, scratch_end);
    if (dec < 0)
      return -1;
    output += decoded_size;
    if (output != output_end)
      return -1;
    return dec;
  }
}

int Krak_DecodeRLE(const byte *src, const size_t src_size, byte *dst, const int dst_size, uint8 *scratch, const uint8 *scratch_end) {
  if (src_size <= 1) {
    if (src_size != 1)
      return -1;
    memset(dst, src[0], dst_size);
    return 1;
  }
  const uint8 *dst_end = dst + dst_size;
  const uint8 *cmd_ptr = src + 1, *cmd_ptr_end = src + src_size;
  // Unpack the first X bytes of the command buffer?
  if (src[0]) {
    uint8 *dst_ptr = scratch;
    int dec_size;
    const int n = Kraken_DecodeBytes(&dst_ptr, src, src + src_size, &dec_size, scratch_end - scratch, true, scratch, scratch_end);
    if (n <= 0)
      return -1;
    const int cmd_len = src_size - n + dec_size;
    if (cmd_len > scratch_end - scratch)
      return -1;
    memcpy(dst_ptr + dec_size, src + n, src_size - n);
    cmd_ptr = dst_ptr;
    cmd_ptr_end = &dst_ptr[cmd_len];
  }

  int rle_byte = 0;

  while (cmd_ptr < cmd_ptr_end) {
    const uint32 cmd = cmd_ptr_end[-1];
    if (cmd - 1 >= 0x2f) {
      cmd_ptr_end--;
      const uint32 bytes_to_copy = -1 - cmd & 0xF;
      const uint32 bytes_to_rle = cmd >> 4;
      if (dst_end - dst < bytes_to_copy + bytes_to_rle || cmd_ptr_end - cmd_ptr < bytes_to_copy)
        return -1;
      memcpy(dst, cmd_ptr, bytes_to_copy);
      cmd_ptr += bytes_to_copy;
      dst += bytes_to_copy;
      memset(dst, rle_byte, bytes_to_rle);
      dst += bytes_to_rle;
    } else if (cmd >= 0x10) {
      const uint32 data = *(uint16*)(cmd_ptr_end - 2) - 4096;
      cmd_ptr_end -= 2;
      const uint32 bytes_to_copy = data & 0x3F;
      const uint32 bytes_to_rle = data >> 6;
      if (dst_end - dst < bytes_to_copy + bytes_to_rle || cmd_ptr_end - cmd_ptr < bytes_to_copy)
        return -1;
      memcpy(dst, cmd_ptr, bytes_to_copy);
      cmd_ptr += bytes_to_copy;
      dst += bytes_to_copy;
      memset(dst, rle_byte, bytes_to_rle);
      dst += bytes_to_rle;
    } else if (cmd == 1) {
      rle_byte = *cmd_ptr++;
      cmd_ptr_end--;
    } else if (cmd >= 9) {
      const uint32 bytes_to_rle = (*(uint16*)(cmd_ptr_end - 2) - 0x8ff) * 128;
      cmd_ptr_end -= 2;
      if (dst_end - dst < bytes_to_rle)
        return -1;
      memset(dst, rle_byte, bytes_to_rle);
      dst += bytes_to_rle;
    } else {
      const uint32 bytes_to_copy = (*(uint16*)(cmd_ptr_end - 2) - 511) * 64;
      cmd_ptr_end -= 2;
      if (cmd_ptr_end - cmd_ptr < bytes_to_copy || dst_end - dst < bytes_to_copy)
        return -1;
      memcpy(dst, cmd_ptr, bytes_to_copy);
      dst += bytes_to_copy;
      cmd_ptr += bytes_to_copy;
    }
  }
  if (cmd_ptr_end != cmd_ptr)
    return -1;

  if (dst != dst_end)
    return -1;

  return src_size;
}

struct TansData {
  uint32 A_used;
  uint32 B_used;
  uint8 A[256];
  uint32 B[256];
};

template<typename T> void SimpleSort(T *p, T *pend) {
  if (p != pend) {
    for (T *lp = p + 1, *rp; lp != pend; lp++) {
      T t = lp[0];
      for (rp = lp; rp > p && t < rp[-1]; rp--)
        rp[0] = rp[-1];
      rp[0] = t;
    }
  }
}

bool Tans_DecodeTable(BitReader *bits, const int L_bits, TansData *tans_data) {
  BitReader_Refill(bits);
  if (BitReader_ReadBitNoRefill(bits)) {
    const int Q = BitReader_ReadBitsNoRefill(bits, 3);
    const int num_symbols = BitReader_ReadBitsNoRefill(bits, 8) + 1;
    if (num_symbols < 2)
      return false;
    int fluff = BitReader_ReadFluff(bits, num_symbols);
    const int total_rice_values = fluff + num_symbols;
    uint8 rice[512 + 16];
    BitReader2 br2;

    // another bit reader...
    br2.p = bits->p - (static_cast<uint>(24 - bits->bitpos + 7) >> 3);
    br2.p_end = bits->p_end;
    br2.bitpos = bits->bitpos - 24 & 7;
    
    if (!DecodeGolombRiceLengths(rice, total_rice_values, &br2))
      return false;
    memset(rice + total_rice_values, 0, 16);

    // Switch back to other bitreader impl
    bits->bitpos = 24;
    bits->p = br2.p;
    bits->bits = 0;
    BitReader_Refill(bits);
    bits->bits <<= br2.bitpos;
    bits->bitpos += br2.bitpos;

    HuffRange range[133];
    fluff = Huff_ConvertToRanges(range, num_symbols, fluff, &rice[num_symbols], bits);
    if (fluff < 0)
      return false;

    BitReader_Refill(bits);

    const uint32 L = 1 << L_bits;
    const uint8 *cur_rice_ptr = rice;
    int average = 6;
    int somesum = 0;
    uint8 *tanstable_A = tans_data->A;
    uint32 *tanstable_B = tans_data->B;

    for (int ri = 0; ri < fluff; ri++) {
      int symbol = range[ri].symbol;
      int num = range[ri].num;
      do {
        BitReader_Refill(bits);

        const int nextra = Q + *cur_rice_ptr++;
        if (nextra > 15)
          return false;
        int v = BitReader_ReadBitsNoRefillZero(bits, nextra) + (1 << nextra) - (1 << Q);

        const int average_div4 = average >> 2;
        int limit = 2 * average_div4;
        if (v <= limit)
          v = average_div4 + (-(v & 1) ^ (static_cast<uint32>(v) >> 1));
        if (limit > v)
          limit = v;  
        v += 1;
        average += limit - average_div4;
        *tanstable_A = symbol;
        *tanstable_B = (symbol << 16) + v;
        tanstable_A += v == 1;
        tanstable_B += v >= 2;
        somesum += v;
        symbol += 1;
      } while (--num);
    }
    tans_data->A_used = tanstable_A - tans_data->A;
    tans_data->B_used = tanstable_B - tans_data->B;
    if (somesum != L)
      return false;

    return true;
  } else {
    bool seen[256];
    memset(seen, 0, sizeof(seen));
    const uint32 L = 1 << L_bits;

    int count = BitReader_ReadBitsNoRefill(bits, 3) + 1;

    const int bits_per_sym = BSR(L_bits) + 1;
    const int max_delta_bits = BitReader_ReadBitsNoRefill(bits, bits_per_sym);

    if (max_delta_bits == 0 || max_delta_bits > L_bits)
      return false;

    uint8 *tanstable_A = tans_data->A;
    uint32 *tanstable_B = tans_data->B;

    int weight = 0;
    int total_weights = 0;

    do {
      BitReader_Refill(bits);

      const int sym = BitReader_ReadBitsNoRefill(bits, 8);
      if (seen[sym])
        return false;

      const int delta = BitReader_ReadBitsNoRefill(bits, max_delta_bits);

      weight += delta;

      if (weight == 0)
        return false;

      seen[sym] = true;
      if (weight == 1) {
        *tanstable_A++ = sym;
      } else {
        *tanstable_B++ = (sym << 16) + weight;
      }

      total_weights += weight;
    } while (--count);

    BitReader_Refill(bits);

    const int sym = BitReader_ReadBitsNoRefill(bits, 8);
    if (seen[sym])
      return false;

    if (L - total_weights < weight || L - total_weights <= 1)
      return false;

    *tanstable_B++ = (sym << 16) + (L - total_weights);

    tans_data->A_used = tanstable_A - tans_data->A;
    tans_data->B_used = tanstable_B - tans_data->B;

    SimpleSort(tans_data->A, tanstable_A);
    SimpleSort(tans_data->B, tanstable_B);
    return true;
  }
}

struct TansLutEnt {
  uint32 x;
  uint8 bits_x;
  uint8 symbol;
  uint16 w;
};

void Tans_InitLut(const TansData *tans_data, const int L_bits, TansLutEnt *lut) {
  TansLutEnt *pointers[4];

  const int L = 1 << L_bits;
  const int a_used = tans_data->A_used;

  const uint slots_left_to_alloc = L - a_used;

  const uint sa = slots_left_to_alloc >> 2;
  pointers[0] = lut;
  uint sb = sa + ((slots_left_to_alloc & 3) > 0);
  pointers[1] = lut + sb;
  sb += sa + ((slots_left_to_alloc & 3) > 1);
  pointers[2] = lut + sb;
  sb += sa + ((slots_left_to_alloc & 3) > 2);
  pointers[3] = lut + sb;

  // Setup the single entrys with weight=1
  {
    TansLutEnt *lut_singles = lut + slots_left_to_alloc, le;
    le.w = 0;
    le.bits_x = L_bits;
    le.x = (1 << L_bits) - 1;
    for (int i = 0; i < a_used; i++) {
      lut_singles[i] = le;
      lut_singles[i].symbol = tans_data->A[i];
    }
  }

  // Setup the entrys with weight >= 2
  int weights_sum = 0;
  for (int i = 0; i < tans_data->B_used; i++) {
    const int weight = tans_data->B[i] & 0xffff;
    const int symbol = tans_data->B[i] >> 16;
    if (weight > 4) {
      const uint32 sym_bits = BSR(weight);
      int Z = L_bits - sym_bits;
      TansLutEnt le;
      le.symbol = symbol;
      le.bits_x = Z;
      le.x = (1 << Z) - 1;
      le.w = L - 1 & (weight << Z);
      int what_to_add = 1 << Z;
      int X = (1 << (sym_bits + 1)) - weight;

      for (int j = 0; j < 4; j++) {
        TansLutEnt *dst = pointers[j];

        const int Y = (weight + (weights_sum - j - 1 & 3)) >> 2;
        if (X >= Y) {
          for(int n = Y; n; n--) {
            *dst++ = le;
            le.w += what_to_add;
          }
          X -= Y;
        } else {
          for (int n = X; n; n--) {
            *dst++ = le;
            le.w += what_to_add;
          }
          Z--;

          what_to_add >>= 1;
          le.bits_x = Z;
          le.w = 0;
          le.x >>= 1;
          for (int n = Y - X; n; n--) {
            *dst++ = le;
            le.w += what_to_add;
          }
          X = weight;
        }
        pointers[j] = dst;
      }
    } else {
      assert(weight > 0);
      uint32 bits = ((1 << weight) - 1) << (weights_sum & 3);
      bits |= bits >> 4;
      int n = weight, ww = weight;
      do {
        const uint32 idx = BSF(bits);
        bits &= bits - 1;
        TansLutEnt *dst = pointers[idx]++;
        dst->symbol = symbol;
        const uint32 weight_bits = BSR(ww);
        dst->bits_x = L_bits - weight_bits;
        dst->x = (1 << (L_bits - weight_bits)) - 1;
        dst->w = L - 1 & (ww++ << (L_bits - weight_bits));
      } while (--n);
    }
    weights_sum += weight;
  }
}

struct TansDecoderParams {
  TansLutEnt *lut;
  uint8 *dst, *dst_end;
  const uint8 *ptr_f, *ptr_b;
  uint32 bits_f, bits_b;
  int bitpos_f, bitpos_b;
  uint32 state_0, state_1, state_2, state_3, state_4;
};

bool Tans_Decode(const TansDecoderParams *params) {
  TansLutEnt *lut = params->lut, *e;
  uint8 *dst = params->dst, *dst_end = params->dst_end;
  const uint8 *ptr_f = params->ptr_f, *ptr_b = params->ptr_b;
  uint32 bits_f = params->bits_f, bits_b = params->bits_b;
  int bitpos_f = params->bitpos_f, bitpos_b = params->bitpos_b;
  uint32 state_0 = params->state_0, state_1 = params->state_1;
  uint32 state_2 = params->state_2, state_3 = params->state_3;
  uint32 state_4 = params->state_4;

  if (ptr_f > ptr_b)
    return false;

#define TANS_FORWARD_BITS()                     \
    bits_f |= *(uint32 *)ptr_f << bitpos_f;     \
    ptr_f += (31 - bitpos_f) >> 3;              \
    bitpos_f |= 24;

#define TANS_FORWARD_ROUND(state)               \
    e = &lut[state];                            \
    *dst++ = e->symbol;                         \
    bitpos_f -= e->bits_x;                      \
    state = (bits_f & e->x) + e->w;             \
    bits_f >>= e->bits_x;                       \
    if (dst >= dst_end)                         \
      break;

#define TANS_BACKWARD_BITS()                    \
    bits_b |= _byteswap_ulong(((uint32 *)ptr_b)[-1]) << bitpos_b;     \
    ptr_b -= (31 - bitpos_b) >> 3;              \
    bitpos_b |= 24;

#define TANS_BACKWARD_ROUND(state)              \
    e = &lut[state];                            \
    *dst++ = e->symbol;                         \
    bitpos_b -= e->bits_x;                      \
    state = (bits_b & e->x) + e->w;             \
    bits_b >>= e->bits_x;                       \
    if (dst >= dst_end)                         \
      break;
  
  if (dst < dst_end) {
    for (;;) {
      TANS_FORWARD_BITS()
      TANS_FORWARD_ROUND(state_0)
      TANS_FORWARD_ROUND(state_1)
      TANS_FORWARD_BITS()
      TANS_FORWARD_ROUND(state_2)
      TANS_FORWARD_ROUND(state_3)
      TANS_FORWARD_BITS()
      TANS_FORWARD_ROUND(state_4)
      TANS_BACKWARD_BITS()
      TANS_BACKWARD_ROUND(state_0)
      TANS_BACKWARD_ROUND(state_1)
      TANS_BACKWARD_BITS()
      TANS_BACKWARD_ROUND(state_2)
      TANS_BACKWARD_ROUND(state_3)
      TANS_BACKWARD_BITS()
      TANS_BACKWARD_ROUND(state_4)
    }
  }

  if (ptr_b - ptr_f + (bitpos_f >> 3) + (bitpos_b >> 3) != 0)
    return false;

  const uint32 states_or = state_0 | state_1 | state_2 | state_3 | state_4;
  if (states_or & ~0xFF)
    return false;

  dst_end[0] = static_cast<uint8>(state_0);
  dst_end[1] = static_cast<uint8>(state_1);
  dst_end[2] = static_cast<uint8>(state_2);
  dst_end[3] = static_cast<uint8>(state_3);
  dst_end[4] = static_cast<uint8>(state_4);
  return true;
}

int Krak_DecodeTans(const byte *src, const size_t src_size, byte *dst, const int dst_size, uint8 *scratch, const uint8 *scratch_end) {
  if (src_size < 8 || dst_size < 5) return -1;
  const uint8 *src_end = src + src_size;
  TansData tans_data;

  BitReader br;
  br.bitpos = 24;
  br.bits = 0;
  br.p = src;
  br.p_end = src_end;
  BitReader_Refill(&br);

  if (BitReader_ReadBitNoRefill(&br)) return -1; // reserved bit
  const int L_bits = BitReader_ReadBitsNoRefill(&br, 2) + 8;
  if (!Tans_DecodeTable(&br, L_bits, &tans_data)) return -1;
  src = br.p - (24 - br.bitpos) / 8;
  if (src >= src_end) return -1;
  const uint32 lut_space_required = (sizeof(TansLutEnt) << L_bits) + 15 &~ 15;
  if (lut_space_required > scratch_end - scratch) return -1;

  TansDecoderParams params;
  params.dst = dst;
  params.dst_end = dst + dst_size - 5;

  params.lut = (TansLutEnt *)ALIGN_POINTER(scratch, 16);
  Tans_InitLut(&tans_data, L_bits, params.lut);

  // Read out the initial state
  const uint32 L_mask = (1 << L_bits) - 1;
  uint32 bits_f = *(uint32*)src;
  src += 4;
  uint32 bits_b = _byteswap_ulong(*(uint32*)(src_end - 4));
  src_end -= 4;
  uint32 bitpos_f = 32, bitpos_b = 32;

  // Read first two.
  params.state_0 = bits_f & L_mask;
  params.state_1 = bits_b & L_mask;
  bits_f >>= L_bits, bitpos_f -= L_bits;
  bits_b >>= L_bits, bitpos_b -= L_bits;

  // Read next two.
  params.state_2 = bits_f & L_mask;
  params.state_3 = bits_b & L_mask;
  bits_f >>= L_bits, bitpos_f -= L_bits;
  bits_b >>= L_bits, bitpos_b -= L_bits;

  // Refill more bits
  bits_f |= *(uint32 *)src << bitpos_f;
  src += (31 - bitpos_f) >> 3;
  bitpos_f |= 24;

  // Read final state variable
  params.state_4 = bits_f & L_mask;
  bits_f >>= L_bits, bitpos_f -= L_bits;

  params.bits_f = bits_f;
  params.ptr_f = src - (bitpos_f >> 3);
  params.bitpos_f = bitpos_f & 7;

  params.bits_b = bits_b;
  params.ptr_b = src_end + (bitpos_b >> 3);
  params.bitpos_b = bitpos_b & 7;

  if (!Tans_Decode(&params))
    return -1;

  return src_size;
}

int Kraken_GetBlockSize(const uint8 *src, const uint8 *src_end, int *dest_size, const int dest_capacity) {
  const byte *src_org = src;
  int src_size, dst_size;

  if (src_end - src < 2)
    return -1; // too few bytes

  const int chunk_type = (src[0] >> 4) & 0x7;
  if (chunk_type == 0) {
    if (src[0] >= 0x80) {
      // In this mode, memcopy stores the length in the bottom 12 bits.
      src_size = ((src[0] << 8) | src[1]) & 0xFFF;
      src += 2;
    } else {
      if (src_end - src < 3)
        return -1; // too few bytes
      src_size = (src[0] << 16) | (src[1] << 8) | src[2];
      if (src_size & ~0x3ffff)
        return -1; // reserved bits must not be set
      src += 3;
    }
    if (src_size > dest_capacity || src_end - src < src_size)
      return -1;
    *dest_size = src_size;
    return src + src_size - src_org;
  }

  if (chunk_type >= 6)
    return -1;

  // In all the other modes, the initial bytes encode
  // the src_size and the dst_size
  if (src[0] >= 0x80) {
    if (src_end - src < 3)
      return -1; // too few bytes

    // short mode, 10 bit sizes
    const uint32 bits = (src[0] << 16) | (src[1] << 8) | src[2];
    src_size = bits & 0x3ff;
    dst_size = src_size + ((bits >> 10) & 0x3ff) + 1;
    src += 3;
  } else {
    // long mode, 18 bit sizes
    if (src_end - src < 5)
      return -1; // too few bytes
    const uint32 bits = (src[1] << 24) | (src[2] << 16) | (src[3] << 8) | src[4];
    src_size = bits & 0x3ffff;
    dst_size = (((bits >> 18) | (src[0] << 14)) & 0x3FFFF) + 1;
    if (src_size >= dst_size)
      return -1;
    src += 5;
  }
  if (src_end - src < src_size || dst_size > dest_capacity)
    return -1;
  *dest_size = dst_size;
  return src_size;
}

int Kraken_DecodeBytes(byte **output, const byte *src, const byte *src_end,
    int *decoded_size, const size_t output_size, const bool force_memmove,
    uint8 *scratch, const uint8 *scratch_end) {
  const byte *src_org = src;
  int src_size, dst_size;

  if (src_end - src < 2) return -1; // too few bytes
  const int chunk_type = src[0] >> 4 & 0x7; //「.???....」
  if (chunk_type == 0) {
    if (src[0] >= 0x80) { //「?.......」
      // In this mode, memcopy stores the length in the bottom 12 bits.
      src_size = (src[0] << 8 | src[1]) & 0xfff; //「1000????|????????」 u12
      src += 2;
    } else {
      if (src_end - src < 3) return -1; // too few bytes
      src_size = src[0] << 16 | src[1] << 8 | src[2]; //「0000..??|????????|????????」 u18
      if (src_size & ~0x3ffff) return -1; // reserved bits must not be set
      src += 3;
    }
    if (src_size > output_size || src_end - src < src_size) return -1;
    *decoded_size = src_size;
    if (force_memmove)
      memmove(*output, src, src_size);
    else
      *output = (byte*)src;
    return src + src_size - src_org;
  }

  // In all the other modes, the initial bytes encode the src_size and the dst_size
  if (src[0] >= 0x80) {
    // short mode, 10 bit sizes
    if (src_end - src < 3) return -1; // too few bytes
    const uint32 bits = src[0] << 16 | src[1] << 8 | src[2];
    src_size = bits & 0x3ff; //「1.......|......??|????????」 u10
    dst_size = src_size + (bits >> 10 & 0x3ff) + 1; //「1...????|??????..|........」 u10
    src += 3;
  } else {
    // long mode, 18 bit sizes
    if (src_end - src < 5) return -1; // too few bytes
    const uint32 bits = src[1] << 24 | src[2] << 16 | src[3] << 8 | src[4];
    src_size = bits & 0x3ffff; //「0.......|........|......??|????????|????????」 u18
    dst_size = ((bits >> 18 | src[0] << 14) & 0x3FFFF) + 1; //「0...????|????????|??????..|........|........」 u18
    if (src_size >= dst_size) return -1;
    src += 5;
  }
  if (src_end - src < src_size || dst_size > output_size) return -1;

  u8 *dst = *output;
  if (dst == scratch) {
    if (scratch_end - scratch < dst_size) return -1;
    scratch += dst_size;
  }

  // printf("%d -> %d (%d)\n", src_size, dst_size, chunk_type);

  int src_used = -1;
  switch (chunk_type) {
  case 2:
  case 4: src_used = Kraken_DecodeBytes_Type12(src, src_size, dst, dst_size, chunk_type >> 1); break;
  case 5: src_used = Krak_DecodeRecursive(src, src_size, dst, dst_size, scratch, scratch_end); break;
  case 3: src_used = Krak_DecodeRLE(src, src_size, dst, dst_size, scratch, scratch_end); break;
  case 1: src_used = Krak_DecodeTans(src, src_size, dst, dst_size, scratch, scratch_end); break;
  }
  if (src_used != src_size) return -1;
  *decoded_size = dst_size;
  return src + src_size - src_org;
}

void CombineScaledOffsetArrays(int *offs_stream, const size_t offs_stream_size, const int scale, const uint8 *low_bits) {
  for (size_t i = 0; i != offs_stream_size; i++)
    offs_stream[i] = scale * offs_stream[i] - low_bits[i];
}

// Unpacks the packed 8 bit offset and lengths into 32 bit.
bool Kraken_UnpackOffsets(const byte *src, const byte *src_end,
                          const byte *packed_offs_stream, const byte *packed_offs_stream_extra, const int packed_offs_stream_size,
                          const int multi_dist_scale,
                          const byte *packed_litlen_stream, const int packed_litlen_stream_size,
                          int *offs_stream, int *len_stream,
                          const bool excess_flag, int excess_bytes) {
  BitReader bits_a, bits_b;
  int i;
  int u32_len_stream_size = 0;
  
  bits_a.bitpos = 24;
  bits_a.bits = 0;
  bits_a.p = src;
  bits_a.p_end = src_end;
  BitReader_Refill(&bits_a);

  bits_b.bitpos = 24;
  bits_b.bits = 0;
  bits_b.p = src_end;
  bits_b.p_end = src;
  BitReader_RefillBackwards(&bits_b);

  if (!excess_flag) {
    if (bits_b.bits < 0x2000) return false;
    int n = 31 - BSR(bits_b.bits);
    bits_b.bitpos += n;
    bits_b.bits <<= n;
    BitReader_RefillBackwards(&bits_b);
    n++;
    u32_len_stream_size = (bits_b.bits >> (32 - n)) - 1;
    bits_b.bitpos += n;
    bits_b.bits <<= n;
    BitReader_RefillBackwards(&bits_b);
  }
  
  if (multi_dist_scale == 0) {
    // Traditional way of coding offsets
    const uint8 *packed_offs_stream_end = packed_offs_stream + packed_offs_stream_size;
    while (packed_offs_stream != packed_offs_stream_end) {
      *offs_stream++ = -static_cast<int32>(BitReader_ReadDistance(&bits_a, *packed_offs_stream++));
      if (packed_offs_stream == packed_offs_stream_end)
        break;
      *offs_stream++ = -static_cast<int32>(BitReader_ReadDistanceB(&bits_b, *packed_offs_stream++));
    }
  } else {
    // New way of coding offsets 
    int *offs_stream_org = offs_stream;
    const uint8 *packed_offs_stream_end = packed_offs_stream + packed_offs_stream_size;
    while (packed_offs_stream != packed_offs_stream_end) {
      uint32 cmd = *packed_offs_stream++;
      if (cmd >> 3 > 26) return false;
      uint32 offs = ((8 + (cmd & 7)) << (cmd >> 3)) | BitReader_ReadMoreThan24Bits(&bits_a, cmd >> 3);
      *offs_stream++ = 8 - static_cast<int32>(offs);
      if (packed_offs_stream == packed_offs_stream_end)
        break;
      cmd = *packed_offs_stream++;
      if (cmd >> 3 > 26) return false;
      offs = ((8 + (cmd & 7)) << (cmd >> 3)) | BitReader_ReadMoreThan24BitsB(&bits_b, cmd >> 3);
      *offs_stream++ = 8 - static_cast<int32>(offs);
    }
    if (multi_dist_scale != 1)
      CombineScaledOffsetArrays(offs_stream_org, offs_stream - offs_stream_org, multi_dist_scale, packed_offs_stream_extra);
  }
  uint32 u32_len_stream_buf[512]; // max count is 128kb / 256 = 512
  if (u32_len_stream_size > 512) return false;
   
  uint32 *u32_len_stream = u32_len_stream_buf,
         *u32_len_stream_end = u32_len_stream_buf + u32_len_stream_size;
  for (i = 0; i + 1 < u32_len_stream_size; i += 2) {
    if (!BitReader_ReadLength(&bits_a, &u32_len_stream[i + 0])) return false;
    if (!BitReader_ReadLengthB(&bits_b, &u32_len_stream[i + 1])) return false;
  }
  if (i < u32_len_stream_size && !BitReader_ReadLength(&bits_a, &u32_len_stream[i + 0])) return false;

  bits_a.p -= (24 - bits_a.bitpos) >> 3;
  bits_b.p += (24 - bits_b.bitpos) >> 3;
  if (bits_a.p != bits_b.p) return false;

  for (i = 0; i < packed_litlen_stream_size; i++) {
    uint32 v = packed_litlen_stream[i];
    if (v == 255)
      v = *u32_len_stream++ + 255;
    len_stream[i] = v + 3;
  }
  if (u32_len_stream != u32_len_stream_end) return false;

  return true;
}

bool Kraken_ReadLzTable(KrakenDecoder* dec, 
    const int mode, const byte *src, const byte *src_end,
    byte *dst, const u32 dst_size, const int offset,
    KrakenLzTable *lztable
) {
  byte *scratch = dec->scratch + sizeof(KrakenLzTable);
  byte *scratch_end = dec->scratch + dec->scratch_size;
  const size_t scratch_usage = Min(Min(3 * dst_size + 0xd020, 0x6C000), scratch_end - dec->scratch);
  if (scratch_usage < sizeof(KrakenLzTable)) return false;
  scratch_end = dec->scratch + scratch_usage;
  
  byte *out;
  int decode_count;
  byte *packed_offs_stream, *packed_len_stream;

  if (mode > 1 || src_end - src < 13) return false;

  if (offset == 0) {
    COPY_64(dst, src)
    dst += 8;
    src += 8;
  }

  if (*src & 0x80) {
    const uint8 flag = *src++;
    if ((flag & 0xc0) != 0x80)
      return false; // reserved flag set

    return false; // excess bytes not supported
  }

  // Disable no copy optimization if source and dest overlap
  const bool force_copy = dst <= src_end && src <= dst + dst_size;

  // Decode lit stream, bounded by dst_size
  out = scratch;
  int n = Kraken_DecodeBytes(&out, src, src_end, &decode_count, Min(scratch_end - scratch, dst_size),
                             force_copy, scratch, scratch_end);
  if (n < 0)
    return false;
  src += n;
  lztable->lit_stream = out;
  lztable->lit_stream_size = decode_count;
  scratch += decode_count;

  // Decode command stream, bounded by dst_size
  out = scratch;
  n = Kraken_DecodeBytes(&out, src, src_end, &decode_count, Min(scratch_end - scratch, dst_size),
    force_copy, scratch, scratch_end);
  if (n < 0)
    return false;
  src += n;
  lztable->cmd_stream = out;
  lztable->cmd_stream_size = decode_count;
  scratch += decode_count;

  // Check if to decode the multistuff crap
  if (src_end - src < 3)
    return false;

  int offs_scaling = 0;
  uint8 *packed_offs_stream_extra = nullptr;

  if (src[0] & 0x80) {
    // uses the mode where distances are coded with 2 tables
    offs_scaling = src[0] - 127;
    src++;

    packed_offs_stream = scratch;
    n = Kraken_DecodeBytes(&packed_offs_stream, src, src_end, &lztable->offs_stream_size,
                           Min(scratch_end - scratch, lztable->cmd_stream_size), false, scratch, scratch_end);
    if (n < 0)
      return false;
    src += n;
    scratch += lztable->offs_stream_size;

    if (offs_scaling != 1) {
      packed_offs_stream_extra = scratch;
      n = Kraken_DecodeBytes(&packed_offs_stream_extra, src, src_end, &decode_count,
                             Min(scratch_end - scratch, lztable->offs_stream_size), false, scratch, scratch_end);
      if (n < 0 || decode_count != lztable->offs_stream_size)
        return false;
      src += n;
      scratch += decode_count;
    }
  } else {
    // Decode packed offset stream, it's bounded by the command length.
    packed_offs_stream = scratch;
    n = Kraken_DecodeBytes(&packed_offs_stream, src, src_end, &lztable->offs_stream_size,
                           Min(scratch_end - scratch, lztable->cmd_stream_size), false, scratch, scratch_end);
    if (n < 0)
      return false;
    src += n;
    scratch += lztable->offs_stream_size;
  }

  // Decode packed litlen stream. It's bounded by 1/4 of dst_size.
  packed_len_stream = scratch;
  n = Kraken_DecodeBytes(&packed_len_stream, src, src_end, &lztable->len_stream_size,
                         Min(scratch_end - scratch, dst_size >> 2), false, scratch, scratch_end);
  if (n < 0)
    return false;
  src += n;
  scratch += lztable->len_stream_size;

  // Reserve memory for final dist stream
  scratch = ALIGN_POINTER(scratch, 16);
  lztable->offs_stream = (int*)scratch;
  scratch += lztable->offs_stream_size * 4;

  // Reserve memory for final len stream
  scratch = ALIGN_POINTER(scratch, 16);
  lztable->len_stream = (int*)scratch;
  scratch += lztable->len_stream_size * 4;

  if (scratch + 64 > scratch_end)
    return false;

  return Kraken_UnpackOffsets(src, src_end, packed_offs_stream, packed_offs_stream_extra,
                              lztable->offs_stream_size, offs_scaling,
                              packed_len_stream, lztable->len_stream_size,
                              lztable->offs_stream, lztable->len_stream, 0, 0);
}

// Note: may access memory out of bounds on invalid input.
bool Kraken_ProcessLzRuns_Type0(const KrakenLzTable *lzt, byte *dst, const byte *dst_end, const byte *dst_start) {
  const byte *cmd_stream = lzt->cmd_stream,
             *cmd_stream_end = cmd_stream + lzt->cmd_stream_size;
  const int *len_stream = lzt->len_stream;
  const int *len_stream_end = lzt->len_stream + lzt->len_stream_size;
  const byte *lit_stream = lzt->lit_stream;
  const byte *lit_stream_end = lzt->lit_stream + lzt->lit_stream_size;
  const int *offs_stream = lzt->offs_stream;
  const int *offs_stream_end = lzt->offs_stream + lzt->offs_stream_size;
  const byte *copyfrom;
  uint32 final_len;
  int32 offset;
  int32 recent_offs[7];
  int32 last_offset;

  recent_offs[3] = -8;
  recent_offs[4] = -8;
  recent_offs[5] = -8;
  last_offset = -8;

  while (cmd_stream < cmd_stream_end) {
    const uint32 f = *cmd_stream++;
    uint32 litlen = f & 3;
    const uint32 offs_index = f >> 6;
    uint32 matchlen = (f >> 2) & 0xF;

    // use cmov
    uint32 next_long_length = *len_stream;
    const int *next_len_stream = len_stream + 1;

    len_stream = litlen == 3 ? next_len_stream : len_stream;
    litlen = litlen == 3 ? next_long_length : litlen;
    recent_offs[6] = *offs_stream;

    COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
    if (litlen > 8) {
      COPY_64_ADD(dst + 8, lit_stream + 8, &dst[last_offset + 8]);
      if (litlen > 16) {
        COPY_64_ADD(dst + 16, lit_stream + 16, &dst[last_offset + 16]);
        if (litlen > 24) {
          do {
            COPY_64_ADD(dst + 24, lit_stream + 24, &dst[last_offset + 24]);
            litlen -= 8;
            dst += 8;
            lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;

    offset = recent_offs[offs_index + 3];
    recent_offs[offs_index + 3] = recent_offs[offs_index + 2];
    recent_offs[offs_index + 2] = recent_offs[offs_index + 1];
    recent_offs[offs_index + 1] = recent_offs[offs_index + 0];
    recent_offs[3] = offset;
    last_offset = offset;

    offs_stream = (int*)((intptr_t)offs_stream + (offs_index + 1 & 4));

    if (static_cast<uintptr_t>(offset) < static_cast<uintptr_t>(dst_start - dst))
      return false; // offset out of bounds

    copyfrom = dst + offset;
    if (matchlen != 15) {
      COPY_64(dst, copyfrom)
      COPY_64(dst + 8, copyfrom + 8)
      dst += matchlen + 2;
    } else {
      matchlen = 14 + *len_stream++; // why is the value not 16 here, the above case copies up to 16 bytes.
      if (static_cast<uintptr_t>(matchlen) >static_cast<uintptr_t>(dst_end - dst))
        return false; // copy length out of bounds
      COPY_64(dst, copyfrom)
      COPY_64(dst + 8, copyfrom + 8)
      COPY_64(dst + 16, copyfrom + 16)
      do {
        COPY_64(dst + 24, copyfrom + 24)
        matchlen -= 8;
        dst += 8;
        copyfrom += 8;
      } while (matchlen > 24);
      dst += matchlen;
    }
  }

  // check for incorrect input
  if (offs_stream != offs_stream_end || len_stream != len_stream_end)
    return false;

  final_len = dst_end - dst;
  if (final_len != lit_stream_end - lit_stream)
    return false;

  if (final_len >= 8) {
    do {
      COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
      dst += 8, lit_stream += 8, final_len -= 8;
    } while (final_len >= 8);
  }
  if (final_len > 0) {
    do {
      *dst = *lit_stream++ + dst[last_offset];
    } while (dst++, --final_len);
  }
  return true;
}

// Note: may access memory out of bounds on invalid input.
bool Kraken_ProcessLzRuns_Type1(const KrakenLzTable *lzt, byte *dst, const byte *dst_end, const byte *dst_start) {
  const byte *cmd_stream = lzt->cmd_stream, 
             *cmd_stream_end = cmd_stream + lzt->cmd_stream_size;
  const int *len_stream = lzt->len_stream;
  const int *len_stream_end = lzt->len_stream + lzt->len_stream_size;
  const byte *lit_stream = lzt->lit_stream;
  const byte *lit_stream_end = lzt->lit_stream + lzt->lit_stream_size;
  const int *offs_stream = lzt->offs_stream;
  const int *offs_stream_end = lzt->offs_stream + lzt->offs_stream_size;
  const byte *copyfrom;
  uint32 final_len;
  int32 offset;
  int32 recent_offs[7];

  recent_offs[3] = -8;
  recent_offs[4] = -8;
  recent_offs[5] = -8;

  while (cmd_stream < cmd_stream_end) {
    const uint32 f = *cmd_stream++;
    uint32 litlen = f & 3;
    const uint32 offs_index = f >> 6;
    uint32 matchlen = (f >> 2) & 0xF;
  
    // use cmov
    uint32 next_long_length = *len_stream;
    const int *next_len_stream = len_stream + 1;

    len_stream = litlen == 3 ? next_len_stream : len_stream; 
    litlen = litlen == 3 ? next_long_length : litlen;
    recent_offs[6] = *offs_stream;

    COPY_64(dst, lit_stream)
    if (litlen > 8) {
      COPY_64(dst + 8, lit_stream + 8)
      if (litlen > 16) {
        COPY_64(dst + 16, lit_stream + 16)
        if (litlen > 24) {
          do {
            COPY_64(dst + 24, lit_stream + 24)
            litlen -= 8;
            dst += 8;
            lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;

    offset = recent_offs[offs_index + 3];
    recent_offs[offs_index + 3] = recent_offs[offs_index + 2];
    recent_offs[offs_index + 2] = recent_offs[offs_index + 1];
    recent_offs[offs_index + 1] = recent_offs[offs_index + 0];
    recent_offs[3] = offset;
    
    offs_stream = (int*)((intptr_t)offs_stream + (offs_index + 1 & 4));

    if (static_cast<uintptr_t>(offset) < static_cast<uintptr_t>(dst_start - dst))
      return false; // offset out of bounds

    copyfrom = dst + offset;
    if (matchlen != 15) {
      COPY_64(dst, copyfrom)
      COPY_64(dst + 8, copyfrom + 8)
      dst += matchlen + 2;
    } else {
      matchlen = 14 + *len_stream++; // why is the value not 16 here, the above case copies up to 16 bytes.
      if (static_cast<uintptr_t>(matchlen) > static_cast<uintptr_t>(dst_end - dst))
        return false; // copy length out of bounds
      COPY_64(dst, copyfrom)
      COPY_64(dst + 8, copyfrom + 8)
      COPY_64(dst + 16, copyfrom + 16)
      do {
        COPY_64(dst + 24, copyfrom + 24)
        matchlen -= 8;
        dst += 8;
        copyfrom += 8;
      } while (matchlen > 24);
      dst += matchlen;
    }
  }

  // check for incorrect input
  if (offs_stream != offs_stream_end || len_stream != len_stream_end)
    return false;

  final_len = dst_end - dst;
  if (final_len != lit_stream_end - lit_stream)
    return false;

  if (final_len >= 64) {
    do {
      COPY_64_BYTES(dst, lit_stream)
      dst += 64, lit_stream += 64, final_len -= 64;
    } while (final_len >= 64);
  }
  if (final_len >= 8) {
    do {
      COPY_64(dst, lit_stream)
      dst += 8, lit_stream += 8, final_len -= 8;
    } while (final_len >= 8);
  }
  if (final_len > 0) {
    do {
      *dst++ = *lit_stream++;
    } while (--final_len);
  }
  return true;
}

bool Kraken_ProcessLzRuns(KrakenDecoder *dec, const int mode, int src_used, const u32 dst_size) {
  KrakenLzTable *lztable = (KrakenLzTable *)dec->scratch;
  if (!Kraken_ReadLzTable(dec, mode, dec->src, dec->src + src_used, dec->dst1, dst_size, dec->offset, lztable)) return false;
  
  const byte *dst_end = dec->dst1 + dst_size;
  if (mode == 1)
    return Kraken_ProcessLzRuns_Type1(lztable, dec->dst1 + (dec->offset == 0 ? 8 : 0), dst_end, dec->dst);
  if (mode == 0)
    return Kraken_ProcessLzRuns_Type0(lztable, dec->dst1 + (dec->offset == 0 ? 8 : 0), dst_end, dec->dst);
  return false;
}

// Decode one 256kb big quantum block. It's divided into two 128k blocks
// internally that are compressed separately but with a shared history.
int Kraken_DecodeQuantum(KrakenDecoder *dec, const uint32 compressed_size) {
  const auto dst_bytes_left = std::min((u32)0x40000, dec->dst_len);
  const byte *dst_end = dec->dst1 + dst_bytes_left;
  const byte *src_end = dec->src + compressed_size;
  byte *scratch_end = dec->scratch + dec->scratch_size;
  const byte *src_in = dec->src;
  int src_used, written_bytes;

  while (dst_end != dec->dst1) {
    const u32 dst_count = std::min((int)(dst_end - dec->dst1), 0x20000);
    if (src_end - dec->src < 4) return -1;
    const int chunk_hdr = dec->src[2] | dec->src[1] << 8 | dec->src[0] << 16;
    if (chunk_hdr & 0x800000) {
      dec->src += 3;
      src_used = chunk_hdr & 0x7ffff;
      const int mode = (chunk_hdr >> 19) & 0xf;
      if (src_end - dec->src < src_used) return -1;
      if (src_used >= dst_count) {
        if (src_used > dst_count || mode != 0) return -1;
        memmove(dec->dst1, dec->src, dst_count);
      } else {
        if (!Kraken_ProcessLzRuns(dec, mode, src_used, dst_count)) return -1;
      }
    } else {
      // Stored as entropy without any match copying.
      src_used = Kraken_DecodeBytes(&dec->dst1, dec->src, src_end, &written_bytes, dst_count, false, dec->scratch, scratch_end);
      if (src_used < 0 || written_bytes != dst_count) return -1;
    }
    dec->src += src_used;
    dec->dst1 += dst_count;
  }
  return dec->src - src_in;
}

namespace {
bool read_bytes(KrakenDecoder *dec, byte **out, u32* decode_count,
    const byte* src_end, const u32 dst_size, byte **scratch0, const byte *scratch_end) {
  const int n = Kraken_DecodeBytes(out, dec->src, src_end,
    (int*)decode_count, dst_size, false, *scratch0, scratch_end);
  if (n < 0) return false;
  dec->src += n;
  return true;
}

struct LeviathanLzTable {
  int *offs_stream;
  u32 offs_stream_size;
  int *len_stream;
  u32 len_stream_size;
  int lit_stream_size[16];
  int lit_stream_total;
  int cmd_stream_size;
  uint8 *lit_stream[16];
  uint8 *multi_cmd_ptr[8];
  uint8 *multi_cmd_end[8];
  uint8 *cmd_stream;
};

bool Leviathan_ReadLzTable(
  KrakenDecoder *dec,
  const int chunk_type,
  const int src_used,
  const int dst_size,
  LeviathanLzTable *lztable
) {
  const size_t scratch_usage = Min(Min(3 * dst_size + 0xd020, 0x6c000), dec->scratch_size);
  if (scratch_usage < sizeof(LeviathanLzTable)) return false;
  byte *scratch_end = dec->scratch + scratch_usage;
  // const byte *src = dec->src;
  const byte *src_end = dec->src + src_used;
  byte *dst = dec->dst1;
  const int offset = dec->offset;
  byte *scratch = dec->scratch + sizeof(LeviathanLzTable);
  byte *packed_offs_stream, *packed_len_stream, *out;
  int decode_count, n;

  if (chunk_type > 5 || src_end - dec->src < 13) return false;
  if (offset == 0) {
    COPY_64(dst, dec->src)
    dst += 8;
    dec->src += 8;
  }

  int offs_scaling = 0;
  uint8 *packed_offs_stream_extra = nullptr;
  const int offs_stream_limit = dst_size / 3;

  u32 output_size = Min(scratch_end - scratch, offs_stream_limit);
  if (!(dec->src[0] & 0x80)) {
    // Decode packed offset stream, it's bounded by the command length.
    read_bytes(dec, &packed_offs_stream, &lztable->offs_stream_size, src_end, output_size, &scratch, scratch_end);
  } else {
    // uses the mode where distances are coded with 2 tables
    // and the transformation offs * scaling + low_bits
    offs_scaling = dec->src[0] - 127;
    dec->src++;

    read_bytes(dec, &packed_offs_stream, &lztable->offs_stream_size, src_end, output_size, &scratch, scratch_end);

    if (offs_scaling != 1) {
      u32 decode_count1;
      read_bytes(dec, &packed_offs_stream_extra, &decode_count1, src_end, output_size, &scratch, scratch_end);
      if (decode_count1 != lztable->offs_stream_size) return false;
    }
  }

  // Decode packed litlen stream. It's bounded by 1/5 of dst_size.
  output_size = Min(scratch_end - scratch, dst_size / 5);
  read_bytes(dec, &packed_len_stream, &lztable->len_stream_size, src_end, output_size, &scratch, scratch_end);
  // packed_len_stream = scratch;
  // n = Kraken_DecodeBytes(&packed_len_stream, dec->src, src_end, &lztable->len_stream_size,
  //                        Min(scratch_end - scratch, dst_size / 5), false, scratch, scratch_end);
  // if (n < 0) return false;
  // dec->src += n;
  // scratch += lztable->len_stream_size;

  // Reserve memory for final dist stream
  scratch = ALIGN_POINTER(scratch, 16);
  lztable->offs_stream = (int*)scratch;
  scratch += lztable->offs_stream_size * 4;

  // Reserve memory for final len stream
  scratch = ALIGN_POINTER(scratch, 16);
  lztable->len_stream = (int*)scratch;
  scratch += lztable->len_stream_size * 4;

  if (scratch > scratch_end) return false;

  if (chunk_type <= 1) {
    // Decode lit stream, bounded by dst_size
    out = scratch;
    n = Kraken_DecodeBytes(&out, dec->src, src_end, &decode_count, Min(scratch_end - scratch, dst_size),
                           true, scratch, scratch_end);
    if (n < 0) return false;
    dec->src += n;
    lztable->lit_stream[0] = out;
    lztable->lit_stream_size[0] = decode_count;
  } else {
    const int array_count = chunk_type == 2 ? 2 :
                      chunk_type == 3 ? 4 : 16;
    n = Kraken_DecodeMultiArray(dec->src, src_end, scratch, scratch_end, lztable->lit_stream,
                                lztable->lit_stream_size, array_count, &decode_count,
                                true, scratch, scratch_end);
    if (n < 0) return false;
    dec->src += n;
  }
  scratch += decode_count;
  lztable->lit_stream_total = decode_count;
  if (dec->src >= src_end) return false;

  if (!(dec->src[0] & 0x80)) {
    // Decode command stream, bounded by dst_size
    out = scratch;
    n = Kraken_DecodeBytes(&out, dec->src, src_end, &decode_count, Min(scratch_end - scratch, dst_size),
                           true, scratch, scratch_end);
    if (n < 0) return false;
    dec->src += n;
    lztable->cmd_stream = out;
    lztable->cmd_stream_size = decode_count;
    scratch += decode_count;
  } else {
    if (dec->src[0] != 0x83) return false;
    dec->src++;
    int multi_cmd_lens[8];
    n = Kraken_DecodeMultiArray(dec->src, src_end, scratch, scratch_end, lztable->multi_cmd_ptr,
                                multi_cmd_lens, 8, &decode_count, true, scratch, scratch_end);
    if (n < 0) return false;
    dec->src += n;
    for (size_t i = 0; i < 8; i++)
      lztable->multi_cmd_end[i] = lztable->multi_cmd_ptr[i] + multi_cmd_lens[i];

    lztable->cmd_stream = nullptr;
    lztable->cmd_stream_size = decode_count;
    scratch += decode_count;
  }
  if (dst_size > scratch_end - scratch) return false;

  return Kraken_UnpackOffsets(dec->src, src_end, packed_offs_stream, packed_offs_stream_extra,
                              lztable->offs_stream_size, offs_scaling,
                              packed_len_stream, lztable->len_stream_size,
                              lztable->offs_stream, lztable->len_stream, 0, 0);
}

#define finline __forceinline

struct LeviathanModeRaw {
  const uint8 *lit_stream;

  finline LeviathanModeRaw(const LeviathanLzTable *lzt, uint8 *dst_start) : lit_stream(lzt->lit_stream[0]) {
  }
  
  finline bool CopyLiterals(const uint32 cmd, uint8 *&dst, const int *&len_stream, const uint8 *match_zone_end, size_t last_offset) {
    uint32 litlen = (cmd >> 3) & 3;
    // use cmov
    const uint32 len_stream_value = *len_stream & 0xffffff;
    const int *next_len_stream = len_stream + 1;
    len_stream = litlen == 3 ? next_len_stream : len_stream;
    litlen = litlen == 3 ? len_stream_value : litlen;
    COPY_64(dst, lit_stream)
    if (litlen > 8) {
      COPY_64(dst + 8, lit_stream + 8)
      if (litlen > 16) {
        COPY_64(dst + 16, lit_stream + 16)
        if (litlen > 24) {
          if (litlen > match_zone_end - dst)
            return false;  // out of bounds
          do {
            COPY_64(dst + 24, lit_stream + 24)
            litlen -= 8, dst += 8, lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;
    return true;
  }

  finline void CopyFinalLiterals(uint32 final_len, uint8 *&dst, size_t last_offset) {
    if (final_len >= 64) {
      do {
        COPY_64_BYTES(dst, lit_stream)
        dst += 64; lit_stream += 64; final_len -= 64;
      } while (final_len >= 64);
    }
    if (final_len >= 8) {
      do {
        COPY_64(dst, lit_stream)
        dst += 8; lit_stream += 8; final_len -= 8;
      } while (final_len >= 8);
    }
    if (final_len > 0) {
      do {
        *dst++ = *lit_stream++;
      } while (--final_len);
    }
  }
};

struct LeviathanModeSub {
  const uint8 *lit_stream;

  finline LeviathanModeSub(const LeviathanLzTable *lzt, uint8 *dst_start) : lit_stream(lzt->lit_stream[0]) {
  }

  finline bool CopyLiterals(const uint32 cmd, uint8 *&dst, const int *&len_stream, const uint8 *match_zone_end, const size_t last_offset) {
    uint32 litlen = (cmd >> 3) & 3;
    // use cmov
    const uint32 len_stream_value = *len_stream & 0xffffff;
    const int *next_len_stream = len_stream + 1;
    len_stream = litlen == 3 ? next_len_stream : len_stream;
    litlen = litlen == 3 ? len_stream_value : litlen;
    COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
    if (litlen > 8) {
      COPY_64_ADD(dst + 8, lit_stream + 8, &dst[last_offset + 8]);
      if (litlen > 16) {
        COPY_64_ADD(dst + 16, lit_stream + 16, &dst[last_offset + 16]);
        if (litlen > 24) {
          if (litlen > match_zone_end - dst)
            return false;  // out of bounds
          do {
            COPY_64_ADD(dst + 24, lit_stream + 24, &dst[last_offset + 24]);
            litlen -= 8, dst += 8, lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;
    return true;
  }

  finline void CopyFinalLiterals(uint32 final_len, uint8 *&dst, const size_t last_offset) {
    if (final_len >= 8) {
      do {
        COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
        dst += 8; lit_stream += 8; final_len -= 8;
      } while (final_len >= 8);
    }
    if (final_len > 0) {
      do {
        *dst = *lit_stream++ + dst[last_offset];
      } while (dst++, --final_len);
    }
  }
};

struct LeviathanModeLamSub {
  const uint8 *lit_stream, *lam_lit_stream;

  finline LeviathanModeLamSub(const LeviathanLzTable *lzt, uint8 *dst_start) 
    : lit_stream(lzt->lit_stream[0]),
      lam_lit_stream(lzt->lit_stream[1]) {
  }

  finline bool CopyLiterals(const uint32 cmd, uint8 *&dst, const int *&len_stream, const uint8 *match_zone_end, const size_t last_offset) {
    const uint32 lit_cmd = cmd & 0x18;
    if (!lit_cmd) return true;

    uint32 litlen = lit_cmd >> 3;
    // use cmov
    const uint32 len_stream_value = *len_stream & 0xffffff;
    const int *next_len_stream = len_stream + 1;
    len_stream = litlen == 3 ? next_len_stream : len_stream;
    litlen = litlen == 3 ? len_stream_value : litlen;
       
    if (litlen-- == 0) return false; // lamsub mode requires one literal

    dst[0] = *lam_lit_stream++ + dst[last_offset]; dst++;

    COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
    if (litlen > 8) {
      COPY_64_ADD(dst + 8, lit_stream + 8, &dst[last_offset + 8]);
      if (litlen > 16) {
        COPY_64_ADD(dst + 16, lit_stream + 16, &dst[last_offset + 16]);
        if (litlen > 24) {
          if (litlen > match_zone_end - dst)
            return false;  // out of bounds
          do {
            COPY_64_ADD(dst + 24, lit_stream + 24, &dst[last_offset + 24]);
            litlen -= 8; dst += 8; lit_stream += 8;
          } while (litlen > 24);
        }
      }
    }
    dst += litlen;
    lit_stream += litlen;
    return true;
  }

  finline void CopyFinalLiterals(uint32 final_len, uint8 *&dst, const size_t last_offset) {
    dst[0] = *lam_lit_stream++ + dst[last_offset]; dst++;
    final_len -= 1;

    if (final_len >= 8) {
      do {
        COPY_64_ADD(dst, lit_stream, &dst[last_offset]);
        dst += 8; lit_stream += 8; final_len -= 8;
      } while (final_len >= 8);
    }
    if (final_len > 0) {
      do {
        *dst = *lit_stream++ + dst[last_offset];
      } while (dst++, --final_len);
    }
  }
};

struct LeviathanModeSubAnd3 {
  enum { NUM = 4, MASK = NUM - 1};
  const uint8 *lit_stream[NUM];

  finline LeviathanModeSubAnd3(const LeviathanLzTable *lzt, uint8 *dst_start) {
    for (size_t i = 0; i != NUM; i++)
      lit_stream[i] = lzt->lit_stream[-(intptr_t)dst_start + i & MASK];
  }
  finline bool CopyLiterals(const uint32 cmd, uint8 *&dst, const int *&len_stream, const uint8 *match_zone_end, const size_t last_offset) {
    const uint32 lit_cmd = cmd & 0x18;

    if (lit_cmd == 0x18) {
      uint32 litlen = *len_stream++ & 0xffffff;
      if (litlen > match_zone_end - dst)
        return false;
      while (litlen) {
        *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
        dst++; litlen--;
      }
    } else if (lit_cmd) {
      *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
      dst++;
      if (lit_cmd == 0x10) {
        *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
        dst++;
      }
    }
    return true;
  }

  finline void CopyFinalLiterals(uint32 final_len, uint8 *&dst, const size_t last_offset) {
    if (final_len > 0) {
      do {
        *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
      } while (dst++, --final_len);
    }
  }
};

struct LeviathanModeSubAndF {
  enum { NUM = 16, MASK = NUM - 1};
  const uint8 *lit_stream[NUM];
  
  finline LeviathanModeSubAndF(const LeviathanLzTable *lzt, uint8 *dst_start) {
    for(size_t i = 0; i != NUM; i++)
      lit_stream[i] = lzt->lit_stream[-(intptr_t)dst_start + i & MASK];
  }
  finline bool CopyLiterals(const uint32 cmd, uint8 *&dst, const int *&len_stream, const uint8 *match_zone_end, const size_t last_offset) {
    const uint32 lit_cmd = cmd & 0x18;

    if (lit_cmd == 0x18) {
      uint32 litlen = *len_stream++ & 0xffffff;
      if (litlen > match_zone_end - dst)
        return false;
      while (litlen) {
        *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
        dst++; litlen--;
      }
    } else if (lit_cmd) {
      *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
      dst++;
      if (lit_cmd == 0x10) {
        *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
        dst++;
      }
    }
    return true;
  }

  finline void CopyFinalLiterals(uint32 final_len, uint8 *&dst, const size_t last_offset) {
    if (final_len > 0) {
      do {
        *dst = *lit_stream[(uintptr_t)dst & MASK]++ + dst[last_offset];
      } while (dst++, --final_len);
    }
  }
};

struct LeviathanModeO1 {
  const uint8 *lit_streams[16];
  uint8 next_lit[16];
  
  finline LeviathanModeO1(const LeviathanLzTable *lzt, uint8 *dst_start) {
    for (size_t i = 0; i != 16; i++) {
      const uint8 *p = lzt->lit_stream[i];
      next_lit[i] = *p;
      lit_streams[i] = p + 1;
    }
  }

  finline bool CopyLiterals(const uint32 cmd, uint8 *&dst, const int *&len_stream, uint8 *match_zone_end, size_t last_offset) {
    const uint32 lit_cmd = cmd & 0x18;

    if (lit_cmd == 0x18) {
      uint32 litlen = *len_stream++;
      if (static_cast<int32>(litlen) <= 0)
        return false;
      uint context = dst[-1];
      do {
        const size_t slot = context >> 4;
        *dst++ = context = next_lit[slot];
        next_lit[slot] = *lit_streams[slot]++;
      } while (--litlen);
    } else if (lit_cmd) {
      // either 1 or 2
      uint context = dst[-1];
      size_t slot = context >> 4;
      *dst++ = context = next_lit[slot];
      next_lit[slot] = *lit_streams[slot]++;
      if (lit_cmd == 0x10) {
        slot = context >> 4;
        *dst++ = context = next_lit[slot];
        next_lit[slot] = *lit_streams[slot]++;
      }
    }
    return true;
  }

  finline void CopyFinalLiterals(uint32 final_len, uint8 *&dst, size_t last_offset) {
    uint context = dst[-1];
    while (final_len) {
      const size_t slot = context >> 4;
      *dst++ = context = next_lit[slot];
      next_lit[slot] = *lit_streams[slot]++;
      final_len--;
    }
  }
};

template<typename Mode>
bool Leviathan_ProcessLz(LeviathanLzTable *lzt, uint8 *dst,
                         uint8 *dst_start, uint8 *dst_end, uint8 *window_base, bool multi_cmd) {
  const uint8 *cmd_stream = lzt->cmd_stream,
              *cmd_stream_end = cmd_stream + lzt->cmd_stream_size;
  const int *len_stream = lzt->len_stream;
  const int *len_stream_end = len_stream + lzt->len_stream_size;
  
  const int *offs_stream = lzt->offs_stream;
  const int *offs_stream_end = offs_stream + lzt->offs_stream_size;
  const byte *copyfrom;
  uint8 *match_zone_end = dst_end - dst_start >= 16 ? dst_end - 16 : dst_start;

  int32 recent_offs[16];
  recent_offs[8] = recent_offs[9] = recent_offs[10] = recent_offs[11] = -8;
  recent_offs[12] = recent_offs[13] = recent_offs[14] = -8;

  size_t offset = -8;

  Mode mode(lzt, dst_start);

  uint32 cmd_stream_left;
  const uint8 *multi_cmd_stream[8], **cmd_stream_ptr;
  if (multi_cmd) {
    for (size_t i = 0; i != 8; i++)
      multi_cmd_stream[i] = lzt->multi_cmd_ptr[i - (uintptr_t)dst_start & 7];
    cmd_stream_left = lzt->cmd_stream_size;
    cmd_stream_ptr = &multi_cmd_stream[(uintptr_t)dst & 7];
    cmd_stream = *cmd_stream_ptr;
  }

  for(;;) {
    uint32 cmd;
    
    if (!multi_cmd) {
      if (cmd_stream >= cmd_stream_end)
        break;
      cmd = *cmd_stream++;
    } else {
      if (cmd_stream_left == 0)
        break;
      cmd_stream_left--;
      cmd = *cmd_stream;
      *cmd_stream_ptr = cmd_stream + 1;
    }

    const uint32 offs_index = cmd >> 5;
    uint32 matchlen = (cmd & 7) + 2;
    recent_offs[15] = *offs_stream;
    if (!mode.CopyLiterals(cmd, dst, len_stream, match_zone_end, offset)) return false;
    offset = recent_offs[static_cast<size_t>(offs_index) + 8];

    // Permute the recent offsets table
    const simde__m128i temp = simde_mm_loadu_si128((const simde__m128i *)&recent_offs[static_cast<size_t>(offs_index) + 4]);
    simde_mm_storeu_si128((simde__m128i *)&recent_offs[static_cast<size_t>(offs_index) + 1], simde_mm_loadu_si128((const simde__m128i *)&recent_offs[offs_index]));
    simde_mm_storeu_si128((simde__m128i *)&recent_offs[static_cast<size_t>(offs_index) + 5], temp);
    recent_offs[8] = static_cast<int32>(offset);
    offs_stream += offs_index == 7;

    if ((uintptr_t)offset < static_cast<uintptr_t>(window_base - dst))
      return false;  // offset out of bounds
    copyfrom = dst + offset;

    if (matchlen == 9) {
      if (len_stream >= len_stream_end) return false;  // len stream empty
      matchlen = *--len_stream_end + 6;
      COPY_64(dst, copyfrom)
      COPY_64(dst + 8, copyfrom + 8)
      uint8 *next_dst = dst + matchlen;
      if (multi_cmd)
        cmd_stream = *(cmd_stream_ptr = &multi_cmd_stream[(uintptr_t)next_dst & 7]);
      if (matchlen > 16) {
        if (matchlen > static_cast<uintptr_t>(dst_end - 8 - dst)) return false;  // no space in buf
        COPY_64(dst + 16, copyfrom + 16)
        do {
          COPY_64(dst + 24, copyfrom + 24)
          matchlen -= 8; dst += 8; copyfrom += 8;
        } while (matchlen > 24);
      }
      dst = next_dst;
    } else {
      COPY_64(dst, copyfrom)
      dst += matchlen;
      if (multi_cmd)
        cmd_stream = *(cmd_stream_ptr = &multi_cmd_stream[(uintptr_t)dst & 7]);
    }
  }

  if (offs_stream != offs_stream_end || len_stream != len_stream_end) return false; // check for incorrect input

  // copy final literals
  if (dst < dst_end) {
    mode.CopyFinalLiterals(dst_end - dst, dst, offset);
  } else if (dst != dst_end) {
    return false;
  }
  return true;
}

bool Leviathan_ProcessLzRuns(KrakenDecoder *dec, const int chunk_type, const int src_used, const int dst_size) {
  LeviathanLzTable *lzt = (LeviathanLzTable *)dec->scratch;
  if (!Leviathan_ReadLzTable(dec, chunk_type, src_used, dst_size, lzt)) return false;

  uint8 *dst_cur = dec->dst1 + (dec->offset == 0 ? 8 : 0);
  uint8 *dst_end = dec->dst1 + dst_size;
  const bool multi_cmd_mode = lzt->cmd_stream == nullptr;
  switch (chunk_type) {
  case 0: return Leviathan_ProcessLz<LeviathanModeSub>(lzt, dst_cur, dec->dst1, dst_end, dec->dst, multi_cmd_mode);
  case 1: return Leviathan_ProcessLz<LeviathanModeRaw>(lzt, dst_cur, dec->dst1, dst_end, dec->dst, multi_cmd_mode);
  case 2: return Leviathan_ProcessLz<LeviathanModeLamSub>(lzt, dst_cur, dec->dst1, dst_end, dec->dst, multi_cmd_mode);
  case 3: return Leviathan_ProcessLz<LeviathanModeSubAnd3>(lzt, dst_cur, dec->dst1, dst_end, dec->dst, multi_cmd_mode);
  case 4: return Leviathan_ProcessLz<LeviathanModeO1>(lzt, dst_cur, dec->dst1, dst_end, dec->dst, multi_cmd_mode);
  case 5: return Leviathan_ProcessLz<LeviathanModeSubAndF>(lzt, dst_cur, dec->dst1, dst_end, dec->dst, multi_cmd_mode);
  default: return false;
  }
}

// Decode one 256kb big quantum block. It's divided into two 128k blocks
// internally that are compressed separately but with a shared history.
bool Leviathan_DecodeQuantum(KrakenDecoder *dec, const u32 compressed_size) {
  const auto dst_bytes_left = std::min<u32>(dec->dst_len, 0x40000);
  const byte *dst_end = dec->dst1 + dst_bytes_left;
  const byte *src_end = dec->src + compressed_size;
  const byte *scratch_end = dec->scratch + dec->scratch_size;

  while (dst_end - dec->dst1 != 0) {
    int dst_count = dst_end - dec->dst1;
    dst_count = std::min(dst_count, 0x20000);
    if (src_end - dec->src < 4) return false;
    const int chunk_hdr = dec->src[2] | dec->src[1] << 8 | dec->src[0] << 16;
    if (chunk_hdr & 0x800000) {
      dec->src += 3;
      const int src_used = chunk_hdr & 0x7ffff;
      const int mode = (chunk_hdr >> 19) & 0xf;
      if (src_end - dec->src < src_used) return false;
      if (src_used >= dst_count) {
        if (src_used > dst_count || mode != 0) return false;
        memmove(dec->dst1, dec->src, dst_count);
      } else {
        if (!Leviathan_ProcessLzRuns(dec, mode, src_used, dst_count)) return false;
      }
    } else { // Stored as entropy without any match copying.
      int written_bytes;
      const int src_used = Kraken_DecodeBytes(&dec->dst1, dec->src, src_end, &written_bytes, dst_count, false, dec->scratch, scratch_end);
      if (src_used < 0 || written_bytes != dst_count) return false;
    }
    dec->dst1 += dst_count;
  }
  return dec->src == src_end;
}
} // leviathan

namespace {
bool Mermaid_ReadLzTable(KrakenDecoder *dec, const byte *src_end, const u32 dst_size, MermaidLzTable *lz) {
  byte *scratch = dec->scratch + sizeof(MermaidLzTable);
  const byte *scratch_end = scratch + std::min<u32>((dst_size << 1) + 0x4020, 0x40000);
  if (dec->offset == 0) {
    COPY_64(dec->dst1, dec->src)
    dec->dst1 += 8; dec->src += 8;
  }

  lz->cmd_stream_2_offs_end = 0; // better safe than sorry
  // Decode lit stream
  u32 output_size = std::min<u32>((u32)(scratch_end - scratch), dst_size);
  if (!read_bytes(dec, &lz->lit_stream, &lz->cmd_stream_2_offs_end, src_end, output_size, &scratch, scratch_end)) return false;
  lz->lit_stream_end = lz->lit_stream + lz->cmd_stream_2_offs_end;
  // Decode flag stream
  output_size = std::min<u32>((u32)(scratch_end - scratch), dst_size);
  if (!read_bytes(dec, &lz->cmd_stream, &lz->cmd_stream_2_offs_end, src_end, output_size, &scratch, scratch_end)) return false;
  lz->cmd_stream_end = lz->cmd_stream + lz->cmd_stream_2_offs_end;
  if (dst_size <= 0x10000) {
    lz->cmd_stream_2_offs = lz->cmd_stream_2_offs_end;
  } else {
    if (src_end - dec->src < 2) return false;
    lz->cmd_stream_2_offs = *dec->src;
    dec->src += 2;
    if (lz->cmd_stream_2_offs > lz->cmd_stream_2_offs_end) return false;
  }

  if (src_end - dec->src < 2) return false;
  const u32 off16_count = *dec->src;
  dec->src += 2;
  if (off16_count == 0xffff) { // off16 is entropy coded
    uint8 *off16_lo, *off16_hi;
    u32 off16_lo_count, off16_hi_count;
    off16_hi = scratch;
    output_size = std::min<u32>((u32)(scratch_end - scratch), dst_size >> 1);
    if (!read_bytes(dec, &off16_hi, &off16_hi_count, src_end, output_size, &scratch, scratch_end)) return false;

    off16_lo = scratch;
    output_size = std::min<u32>((u32)(scratch_end - scratch), dst_size >> 1);
    if (!read_bytes(dec, &off16_lo, &off16_lo_count, src_end, output_size, &scratch, scratch_end)) return false;

    if (off16_lo_count != off16_hi_count) return false;
    scratch = ALIGN_POINTER(scratch, 2);
    lz->off16_stream = (uint16*)scratch;
    if (scratch + (off16_lo_count << 1) > scratch_end) return false;
    scratch += off16_lo_count << 1;
    lz->off16_stream_end = (uint16*)scratch;

    for (size_t i = 0; i != off16_lo_count; i++) // combine u16 offsets
      lz->off16_stream[i] = off16_lo[i] + off16_hi[i] * 256;
  } else {
    lz->off16_stream = (uint16*)*dec->src;
    dec->src += off16_count << 1;
    lz->off16_stream_end = (uint16*)dec->src;
  }

  if (src_end - dec->src < 3) return false;
  const uint32 tmp = dec->src[0] | dec->src[1] << 8 | dec->src[2] << 16;
  dec->src += 3;
  if (tmp != 0) {
    lz->off32_size_1 = tmp >> 12;
    lz->off32_size_2 = tmp & 0xfff;
    if (lz->off32_size_1 == 0xfff) {
      if (src_end - dec->src < 2) return false;
      lz->off32_size_1 = *(uint16*)dec->src;
      dec->src += 2;
    }
    if (lz->off32_size_2 == 0xfff) {
      if (src_end - dec->src < 2) return false;
      lz->off32_size_2 = *(uint16*)dec->src;
      dec->src += 2;
    }
    if (scratch + ((lz->off32_size_2 + lz->off32_size_1) << 2) + 64 > scratch_end) return false;
    scratch = ALIGN_POINTER(scratch, 4);

    lz->off32_stream_1 = (uint32*)scratch;
    scratch += lz->off32_size_1 << 2;
    // store dummy bytes after for prefetcher.
    ((uint64*)scratch)[0] = ((uint64*)scratch)[1] = ((uint64*)scratch)[2] = ((uint64*)scratch)[3] = 0;
    scratch += 32;

    lz->off32_stream_2 = (uint32*)scratch;
    scratch += lz->off32_size_2 << 2;
    // store dummy bytes after for prefetcher.
    ((uint64*)scratch)[0] = ((uint64*)scratch)[1] = ((uint64*)scratch)[2] = ((uint64*)scratch)[3] = 0;
    scratch += 32;

    // decode far offsets
    for (size_t i = 0; i != lz->off32_size_1; i++) {
      if (src_end - dec->src < 3) return false;
      lz->off32_stream_1[i] = dec->src[0] | dec->src[1] << 8 | dec->src[2] << 16; dec->src += 3;
      if (dec->offset < 0xc00000 - 1 && lz->off32_stream_1[i] >= 0xc00000) {
        if (dec->src == src_end) return false;
        lz->off32_stream_1[i] += *dec->src++ << 22;
      }
      if (lz->off32_stream_1[i] > dec->offset) return false;
    }
    for (size_t i = 0; i != lz->off32_size_2; i++) {
      if (src_end - dec->src < 3) return false;
      lz->off32_stream_2[i] = dec->src[0] | dec->src[1] << 8 | dec->src[2] << 16; dec->src += 3;
      if (dec->offset + 0x10000 < 0xc00000 - 1 && lz->off32_stream_2[i] >= 0xc00000) {
        if (dec->src == src_end) return false;
        lz->off32_stream_2[i] += *dec->src++ << 22;
      }
      if (lz->off32_stream_2[i] > dec->offset + 0x10000) return false;
    }
  } else {
    if (scratch_end - scratch < 32) return false;
    lz->off32_size_1 = lz->off32_size_2 = 0;
    lz->off32_stream_1 = lz->off32_stream_2 = (uint32*)scratch;
    // store dummy bytes after for prefetcher.
    ((uint64*)scratch)[0] = ((uint64*)scratch)[1] = ((uint64*)scratch)[2] = ((uint64*)scratch)[3] = 0;
  }
  lz->length_stream = dec->src;
  return true;
}

const byte *Mermaid_Mode0(const bool mode0, byte *dst, const size_t dst_size,
                          const byte *src_end, MermaidLzTable *lz, int32 *saved_dist, const size_t startoff) {
  const byte *dst_end = dst + dst_size;
  const byte *cmd_stream = lz->cmd_stream;
  const byte *cmd_stream_end = lz->cmd_stream_end;
  const byte *length_stream = lz->length_stream;
  byte *lit_stream = lz->lit_stream;
  const byte *lit_stream_end = lz->lit_stream_end;
  uint16 *off16_stream = lz->off16_stream;
  const uint16 *off16_stream_end = lz->off16_stream_end;
  const uint32 *off32_stream = lz->off32_stream;
  const uint32 *off32_stream_end = lz->off32_stream_end;
  intptr_t recent_offs = *saved_dist;
  const byte *match;
  intptr_t length;
  const byte *dst_begin = dst;

  dst += startoff;

  while (cmd_stream < cmd_stream_end) {
    const uintptr_t cmd = *cmd_stream++;
    if (cmd >= 24) {
      const intptr_t new_dist = *off16_stream;
      const uintptr_t use_distance = (cmd >> 7) - 1;
      const uintptr_t litlen = cmd & 7;
      if (mode0)
        COPY_64_ADD(dst, lit_stream, &dst[recent_offs]);
      else
        COPY_64(dst, lit_stream)
      dst += litlen;
      lit_stream += litlen;
      recent_offs ^= use_distance & (recent_offs ^ -new_dist);
      off16_stream = (uint16*)((uintptr_t)off16_stream + (use_distance & 2));
      match = dst + recent_offs;
      COPY_64(dst, match)
      COPY_64(dst + 8, match + 8)
      dst += (cmd >> 3) & 0xF;
    } else if (cmd > 2) {
      length = cmd + 5;

      if (off32_stream == off32_stream_end) return nullptr;
      match = dst_begin - *off32_stream++;
      recent_offs = match - dst;

      if (dst_end - dst < length) return nullptr;
      COPY_64(dst, match)
      COPY_64(dst + 8, match + 8)
      COPY_64(dst + 16, match + 16)
      COPY_64(dst + 24, match + 24)
      dst += length;
      simde_mm_prefetch((char*)dst_begin - off32_stream[3], SIMDE_MM_HINT_T0);
    } else if (cmd == 0) {
      if (src_end - length_stream == 0) return nullptr;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3) return nullptr;
        length += static_cast<size_t>(*(uint16 *)(length_stream + 1)) * 4;
        length_stream += 2;
      }
      length_stream += 1;

      length += 64;
      if (dst_end - dst < length || lit_stream_end - lit_stream < length) return nullptr;

      do {
        if (mode0) {
          COPY_64_ADD(dst, lit_stream, &dst[recent_offs]);
          COPY_64_ADD(dst + 8, lit_stream + 8, &dst[recent_offs + 8]);
        } else {
          COPY_64(dst, lit_stream)
          COPY_64(dst + 8, lit_stream + 8)
        }
        dst += 16;
        lit_stream += 16;
        length -= 16;
      } while (length > 0);
      dst += length;
      lit_stream += length;
    } else if (cmd == 1) {
      if (src_end - length_stream == 0) return nullptr;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3) return nullptr;
        length += static_cast<size_t>(*(uint16 *)(length_stream + 1)) * 4;
        length_stream += 2;
      }
      length_stream += 1;
      length += 91;

      if (off16_stream == off16_stream_end) return nullptr;
      match = dst - *off16_stream++;
      recent_offs = match - dst;
      do {
        COPY_64(dst, match)
        COPY_64(dst + 8, match + 8)
        dst += 16;
        match += 16;
        length -= 16;
      } while (length > 0);
      dst += length;
    } else /* flag == 2 */ {
      if (src_end - length_stream == 0) return nullptr;
      length = *length_stream;
      if (length > 251) {
        if (src_end - length_stream < 3) return nullptr;
        length += static_cast<size_t>(*(uint16 *)(length_stream + 1)) * 4;
        length_stream += 2;
      }
      length_stream += 1;
      length += 29;
      if (off32_stream == off32_stream_end) return nullptr;
      match = dst_begin - *off32_stream++;
      recent_offs = match - dst;
      do {
        COPY_64(dst, match)
        COPY_64(dst + 8, match + 8)
        dst += 16; match += 16; length -= 16;
      } while (length > 0);
      dst += length;
      simde_mm_prefetch((char*)dst_begin - off32_stream[3], SIMDE_MM_HINT_T0);
    }
  }

  length = dst_end - dst;
  if (length >= 8) {
    do {
      if (mode0)
        COPY_64_ADD(dst, lit_stream, &dst[recent_offs]);
      else
        COPY_64(dst, lit_stream)
      dst += 8; lit_stream += 8; length -= 8;
    } while (length >= 8);
  }
  if (length > 0) {
    do {
      if (mode0) {
        *dst = *lit_stream++ + dst[recent_offs];
        dst++;
      } else {
        *dst++ = *lit_stream++;
      }
    } while (--length);
  }

  *saved_dist = static_cast<int32>(recent_offs);
  lz->length_stream = length_stream;
  lz->off16_stream = off16_stream;
  lz->lit_stream = lit_stream;
  return length_stream;
}

bool Mermaid_ProcessLzRuns(KrakenDecoder *dec, const int mode, const int src_used, u32 dst_size) {
  byte *dst = dec->dst1;
  const u64 offset = dec->offset;
  const auto lz = reinterpret_cast<MermaidLzTable *>(dec->scratch);
  const byte *src_end = dec->src + src_used;
  if (!Mermaid_ReadLzTable(dec, src_end, dst_size, lz)) return false;

  int32 saved_dist = -8;

  // iteration 1
  const u32 dst_size_cur0 = std::min<u32>(dst_size, 0x10000);
  lz->off32_stream = lz->off32_stream_1;
  lz->off32_stream_end = lz->off32_stream_1 + lz->off32_size_1 * 4;
  lz->cmd_stream_end = lz->cmd_stream + lz->cmd_stream_2_offs;
  const byte* src_cur = Mermaid_Mode0(mode == 0, dst, dst_size_cur0, src_end, lz, &saved_dist, offset == 0 ? 8 : 0);
  if (src_cur == nullptr) return false;
  dst_size -= dst_size_cur0;
  if (dst_size == 0) return src_cur == src_end;
  
  // iteration 2
  const u32 dst_size_cur1 = std::min<u32>(dst_size, 0x10000);
  lz->off32_stream = lz->off32_stream_2;
  lz->off32_stream_end = lz->off32_stream_2 + lz->off32_size_2 * 4;
  lz->cmd_stream_end = lz->cmd_stream + lz->cmd_stream_2_offs_end;
  lz->cmd_stream += lz->cmd_stream_2_offs;
  src_cur = Mermaid_Mode0(mode == 0, dst + dst_size_cur0, dst_size_cur1, src_end, lz, &saved_dist, 0);
  if (src_cur == nullptr) return false;
  return src_cur == src_end;
}

bool Mermaid_DecodeQuantum(KrakenDecoder *dec, const u32 compressed_size) {
  const auto dst_bytes_left = std::min<u32>(dec->dst_len, 0x40000);
  const byte *dst_end = dec->dst1 + dst_bytes_left;
  const byte *src_end = dec->src + compressed_size;

  while (dst_end - dec->dst1 != 0) {
    const u32 dst_count = std::min<u32>((u32)(dst_end - dec->dst1), 0x20000);
    if (src_end - dec->src < 4) return false;
    const int chunk_hdr = dec->src[0] << 16 | dec->src[1] << 8 | dec->src[2];
    if (chunk_hdr & 0x800000) { //「?.......|........|........」
      dec->src += 3;
      const int src_used = chunk_hdr & 0x7ffff; //「1....???|????????|????????」 u19
      const int mode = chunk_hdr >> 19 & 0xf; //「1????...|........|........」 u4
      // printf("quantum mode=%d size=%d pos=(%p)\n", mode, src_used, dec->src);
      if (src_end - dec->src < src_used) return false;
      if (src_used >= dst_count) {
        if (src_used > dst_count || mode != 0) return false;
        memmove(dec->dst1, dec->src, dst_count);
      } else { // Tans Lut may need upwards of 16k of temp storage
        if (mode > 1 || src_used < 10) return false;
        if (!Mermaid_ProcessLzRuns(dec, mode, src_used, dst_count)) return false;
      }
    } else { // Stored without any match copying.
      int written_bytes;
      const int src_used = Kraken_DecodeBytes(&dec->dst1, dec->src, src_end, &written_bytes, dst_count, false, dec->scratch, dec->scratch + dec->scratch_size);
      if (src_used < 0 || written_bytes != dst_count) return false;
    }
    dec->dst1 += dst_count;
  }
  return dec->src == src_end;
}
} // mermaid

namespace {
bool step(KrakenDecoder *dec) {
  if ((dec->offset & 0x3ffff) == 0) {
    const int b0 = dec->src[0], b1 = dec->src[1]; dec->src += 2;
    if ((b0 & 0x3f) != 0x0c) return false; // 「..001100|........」
    dec->hdr.uncompressed = (b0 >> 6) & 1; // 「.?......|........」
    dec->hdr.restart_decoder = (b0 >> 7) & 1; // 「?.......|........」
    dec->hdr.decoder_type = b1 & 0x7f; // 「........|.???????」
    dec->hdr.use_checksums = !!(b1 >> 7); // 「........|?.......」
    // printf("hdr: %d %d\n", b0, b1);
  }

  const auto dst_bytes_left = std::min((u32)0x40000, dec->dst_len);

  if (dec->hdr.uncompressed) {
    if (dec->src_len < dst_bytes_left) return false;
    memmove(dec->dst1, dec->src, dst_bytes_left);
    return step_ret(dec, dst_bytes_left, dst_bytes_left);
  }

  const u32 v = dec->src[0] << 16 | dec->src[1] << 8 | dec->src[2]; dec->src += 3;
  const u32 compressed_size = (v & 0x3ffff) + 1;  // 「......??|????????|????????」
  if (compressed_size == 0x40000) {
    printf("memset compressed_size=%u type=%d\n", compressed_size, dec->hdr.decoder_type);
    if (v >> 18 != 1) return false; // 「??????11|11111111|11111111」
    memset(dec->dst1, dec->src[0], dst_bytes_left);
    dec->src += 1;
    return step_ret(dec, 0, dst_bytes_left);
  }
  if (dec->hdr.use_checksums) dec->src += 3; // checksum not implemented
  if (compressed_size > dec->src_len || compressed_size > dst_bytes_left) return false; // Too few bytes in buffer to make any progress?

  int n;
  switch (dec->hdr.decoder_type) {
  case 6: n = Kraken_DecodeQuantum(dec, compressed_size); break;
  case 10: return Mermaid_DecodeQuantum(dec, compressed_size);
  case 12: return Leviathan_DecodeQuantum(dec, compressed_size);
  default: return false;
  }
  if (n != compressed_size) return false;
  return step_ret(dec, n, dst_bytes_left);
}

int kraken_decompress(uint8_t const* src, const size_t src_len, byte *dst, const size_t dst_len) {
  KrakenDecoder *dec = Kraken_Create(src, (u32)src_len, dst, (u32)dst_len);
  while (dec->dst_len != 0)
    if (!step(dec)) goto FAIL;
  if (dec->src_len != 0)
    goto FAIL;
  Kraken_Destroy(dec);
  return (int)(dec->dst1 - dst);
FAIL:
  Kraken_Destroy(dec);
  return -1;
}
}

extern "C" {
    OOZ_DLL_PUBLIC int Ooz_Decompress(uint8_t const* src_buf, const int src_len, uint8_t* dst, const size_t dst_size,
        int, int, int, uint8_t*, size_t, void*, void*, void*, size_t, int) {
        return kraken_decompress(src_buf, src_len, dst, dst_size);
    }
}

#if !OOZ_BUILD_DLL

// The decompressor will write outside the target buffer.
#define SAFE_SPACE 64

void error(const char *s, const char *curfile = NULL) {
  if (curfile)
    fprintf(stderr, "%s: ", curfile);
  fprintf(stderr, "%s\n", s);
  exit(1);
}


byte *load_file(const char *filename, int *size) {
  FILE *f = fopen(filename, "rb");
  if (!f) error("file open error", filename);
  fseek(f, 0, SEEK_END);
  int packed_size = ftell(f);
  fseek(f, 0, SEEK_SET);
  byte *input = new byte[packed_size];
  if (!input) error("memory error", filename);
  if (fread(input, 1, packed_size, f) != packed_size) error("error reading", filename);
  fclose(f);
  *size = packed_size;
  return input;
}

enum {
  kCompressor_Kraken = 8,
  kCompressor_Mermaid = 9,
  kCompressor_Selkie = 11,
  kCompressor_Hydra = 12,
  kCompressor_Leviathan = 13,
};

bool arg_stdout, arg_force, arg_quiet, arg_dll;
int arg_compressor = kCompressor_Kraken, arg_level = 4;
char arg_direction;
const char *verifyfolder;

int ParseCmdLine(int argc, char *argv[]) {
  int i;
  // parse command line
  for (i = 1; i < argc; i++) {
    const char *s = argv[i];
    char c;
    if (*s != '-')
      break;
    if (*++s == '-') {
      if (*++s == 0) {
        i++;
        break;  // --
      }
      // long opts
      if (!strcmp(s, "stdout")) s = "c";
      else if (!strcmp(s, "decompress")) s = "d";
      else if (!strcmp(s, "compress")) s = "z";
      else if (!strncmp(s, "verify=",7)) {
        verifyfolder = s + 7;
        continue;
      } else if (!strcmp(s, "verify")) {
        arg_direction = 't';
        continue;
      } else if (!strcmp(s, "dll")) {
        arg_dll = true;
        continue;
      } else if (!strcmp(s, "kraken")) s = "mk";
      else if (!strcmp(s, "mermaid")) s = "mm";
      else if (!strcmp(s, "selkie")) s = "ms";
      else if (!strcmp(s, "leviathan")) s = "ml";
      else if (!strcmp(s, "hydra")) s = "mh";
      else if (!strncmp(s, "level=", 6)) {
        arg_level = atoi(s + 6);
        continue;
      } else {
        return -1;
      }
    }
    // short opt
    do {
      switch (c = *s++) {
      case 'z':
      case 'd':
      case 'b':
        if (arg_direction)
          return -1;
        arg_direction = c;
        break;
      case 'c':
        arg_stdout = true;
        break;
      case 'f':
        arg_force = true;
        break;
      case 'q':
        arg_quiet = true;
        break;
      case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
        arg_level = c - '0';
        break;
      case 'm':
        c = *s++;
        arg_compressor = (c == 'k') ? kCompressor_Kraken :
                         (c == 'm') ? kCompressor_Mermaid : 
                         (c == 's') ? kCompressor_Selkie : 
                         (c == 'l') ? kCompressor_Leviathan : 
                         (c == 'h') ? kCompressor_Hydra : -1;
        if (arg_compressor < 0)
          return -1;
        break;
      default:
        return -1;
      }
    } while (*s);
  }
  return i;
}

bool Verify(const char *filename, uint8 *output, int outbytes, const char *curfile) {
  int test_size;
  byte *test = load_file(filename, &test_size);
  if (!test) {
    fprintf(stderr, "file open error: %s\n", filename);
    return false;
  }
  if (test_size != outbytes) {
    fprintf(stderr, "%s: ERROR: File size difference: %d vs %d\n", filename, outbytes, test_size);
    return false;
  }
  for (int i = 0; i != test_size; i++) {
    if (test[i] != output[i]) {
      fprintf(stderr, "%s: ERROR: File difference at 0x%x. Was %d instead of %d\n", curfile, i, output[i], test[i]);
      return false;
    }
  }
  return true;
}

#ifndef _MSC_VER
typedef uint64_t LARGE_INTEGER;
void QueryPerformanceCounter(LARGE_INTEGER *a) {
  *a = 0;
}
void QueryPerformanceFrequency(LARGE_INTEGER *a) {
  *a = 1;
}
#define WINAPI
typedef void* HINSTANCE;
HINSTANCE LoadLibraryA(const char *s) { return 0; }
void *GetProcAddress(HINSTANCE h, const char *s) { return 0; }
#endif


typedef int WINAPI OodLZ_CompressFunc(
  int codec, uint8 *src_buf, size_t src_len, uint8 *dst_buf, int level,
  void *opts, size_t offs, size_t unused, void *scratch, size_t scratch_size);
typedef int WINAPI OodLZ_DecompressFunc(uint8 *src_buf, int src_len, uint8 *dst, size_t dst_size,
                                          int fuzz, int crc, int verbose,
                                          uint8 *dst_base, size_t e, void *cb, void *cb_ctx, void *scratch, size_t scratch_size, int threadPhase);

OodLZ_CompressFunc *OodLZ_Compress;
OodLZ_DecompressFunc *OodLZ_Decompress;

void LoadLib() {

#if defined(_M_X64)
#define LIBNAME "oo2core_7_win64.dll"
  char COMPFUNCNAME[] = "XXdleLZ_Compress";
  char DECFUNCNAME[] = "XXdleLZ_Decompress";
  COMPFUNCNAME[0] = DECFUNCNAME[0] = 'O';
  COMPFUNCNAME[1] = DECFUNCNAME[1] = 'o';
#else
#define LIBNAME "oo2core_7_win32.dll"
  char COMPFUNCNAME[] = "_XXdleLZ_Compress@40";
  char DECFUNCNAME[] = "_XXdleLZ_Decompress@56";
  COMPFUNCNAME[1] = DECFUNCNAME[1] = 'O';
  COMPFUNCNAME[2] = DECFUNCNAME[2] = 'o';
#endif
  HINSTANCE mod = LoadLibraryA(LIBNAME);
  OodLZ_Compress = (OodLZ_CompressFunc*)GetProcAddress(mod, COMPFUNCNAME);
  OodLZ_Decompress = (OodLZ_DecompressFunc*)GetProcAddress(mod, DECFUNCNAME);
  if (!OodLZ_Compress || !OodLZ_Decompress)
    error("error loading", LIBNAME);
}

struct CompressOptions;
struct LRMCascade;

int CompressBlock(int codec_id, uint8 *src_in, uint8 *dst_in, int src_size, int level,
                  const CompressOptions *compressopts, uint8 *src_window_base, LRMCascade *lrm);

int main(int argc, char *argv[]) {
  int64_t start, end, freq;
  int argi;

  if (argc < 2 || 
      (argi = ParseCmdLine(argc, argv)) < 0 || 
      argi >= argc ||  // no files
      (arg_direction != 'b' && (argc - argi) > 2) ||  // too many files
      (arg_direction == 't' && (argc - argi) != 2)     // missing argument for verify
      ) {
    fprintf(stderr, "ooz v7.1 - compressor by Rarten\n\n"
      "Usage: ooz [options] input [output]\n"
      " -c --stdout              write to stdout\n"
      " -d --decompress          decompress (default)\n"
      " -z --compress            compress\n"
      " -b                       just benchmark, don't overwrite anything\n"
      " -f                       force overwrite existing file\n"
      " --dll                    decompress with the dll\n"
      " --verify                 decompress and verify that it matches output\n"
      " --verify=<folder>        verify with files in this folder\n"
      " -<1-9> --level=<-4..10>  compression level\n"
      " -m<k>                    [k|m|s|l|h] compressor selection\n"
      " --kraken --mermaid --selkie --leviathan --hydra    compressor selection\n\n"
      "(Warning! not fuzz safe, so please trust the input)\n"
      );
    return 1;
  }
  bool write_mode = (argi + 1 < argc) && (arg_direction != 't' && arg_direction != 'b');

  if (!arg_force && write_mode) {
    struct stat sb;
    if (stat(argv[argi + 1], &sb) >= 0) {
      fprintf(stderr, "file %s already exists, skipping.\n", argv[argi + 1]);
      return 1;
    }
  }

  int nverify = 0;

  for (; argi < argc; argi++) {
    const char *curfile = argv[argi];

    int input_size;
    byte *input = load_file(curfile, &input_size);

    byte *output = NULL;
    int outbytes = 0;

    if (arg_direction == 'z') {
      // compress using the dll
      if (arg_dll)
        LoadLib();
      output = new byte[input_size + 65536];
      if (!output) error("memory error", curfile);
      *(uint64*)output = input_size;
      QueryPerformanceCounter((LARGE_INTEGER*)&start);
      if (arg_dll) {
        outbytes = OodLZ_Compress(arg_compressor, input, input_size, output + 8, arg_level, 0, 0, 0, 0, 0);
      } else {
        outbytes = CompressBlock(arg_compressor, input, output + 8, input_size, arg_level, 0, 0, 0);
      }
      if (outbytes < 0) error("compress failed", curfile);
      outbytes += 8;
      QueryPerformanceCounter((LARGE_INTEGER*)&end);
      QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
      double seconds = (double)(end - start) / freq;
      if (!arg_quiet)
        fprintf(stderr, "%-20s: %8d => %8d (%.2f seconds, %.2f MB/s)\n", argv[argi], input_size, outbytes, seconds, input_size * 1e-6 / seconds);
    } else {
      if (arg_dll)
        LoadLib();

      // stupidly attempt to autodetect if file uses 4-byte or 8-byte header,
      // the previous version of this tool wrote a 4-byte header.
      int hdrsize = *(uint64*)input >= 0x10000000000 ? 4 : 8;
      
      uint64 unpacked_size = (hdrsize == 8) ? *(uint64*)input : *(uint32*)input;
      if (unpacked_size > (hdrsize == 4 ? 52*1024*1024 : 1024 * 1024 * 1024)) 
        error("file too large", curfile);
      output = new byte[unpacked_size + SAFE_SPACE];
      if (!output) error("memory error", curfile);

      QueryPerformanceCounter((LARGE_INTEGER*)&start);

      if (arg_dll) {
        outbytes = OodLZ_Decompress(input + hdrsize, input_size - hdrsize, output, unpacked_size, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      } else {
        outbytes = kraken_decompress(input + hdrsize, input_size - hdrsize, output, unpacked_size);
      }
      if (outbytes != unpacked_size)
        error("decompress error", curfile);
      QueryPerformanceCounter((LARGE_INTEGER*)&end);
      QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
      double seconds = (double)(end - start) / freq;
      if (!arg_quiet)
        fprintf(stderr, "%-20s: %8d => %8d (%.2f seconds, %.2f MB/s)\n", argv[argi], input_size, (int)unpacked_size, seconds, unpacked_size * 1e-6 / seconds);
    }

    if (verifyfolder) {
      // Verify against the file in verifyfolder with the same basename excluding extension
      char buf[1024];
      const char *basename = curfile;
      for(const char *s = curfile; *s; s++)
        if (*s == '/' || *s == '\\')
          basename = s + 1;
      const char *ext = strrchr(basename, '.');
      snprintf(buf, sizeof(buf), "%s/%.*s", verifyfolder, (int)(ext ? (ext - basename) : strlen(basename)), basename);
      if (!Verify(buf, output, outbytes, curfile))
        return 1;
      nverify++;
    }

    if (arg_stdout && arg_direction != 't' && arg_direction != 'b')
      fwrite(output, 1, outbytes, stdout);

    if (write_mode) {
      if (arg_direction == 't') {
        if (!Verify(argv[argi + 1], output, outbytes, curfile))
          return 1;
        fprintf(stderr, "%s: Verify OK\n", curfile);
      } else {
        FILE *f = fopen(argv[argi + 1], "wb");
        if (!f) error("file open for write error", argv[argi + 1]);
        if (fwrite(output, 1, outbytes, f) != outbytes)
          error("file write error", argv[argi + 1]);
        fclose(f);
      }
      break;
    }
    delete[] input;
    delete[] output;
  }

  if (nverify)
    fprintf(stderr, "%d files verified OK!\n", nverify);
  return 0;
}

#endif
