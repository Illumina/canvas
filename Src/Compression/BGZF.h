#pragma once

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

// define the global constants
const int GZIP_WINDOW_BITS = -15;
const int GZIP_ID1   = 31;
const int GZIP_ID2   = 139;
const int CM_DEFLATE = 8;
const int FLG_FEXTRA = 4;
const int OS_UNKNOWN = 255;
const int Z_DEFAULT_MEM_LEVEL = 8;

const int MAX_BLOCK_SIZE = 65536;

const int BGZF_XLEN = 6;
const int BGZF_ID1  = 66;
const int BGZF_ID2  = 67;
const int BGZF_LEN  = 2;

const int BLOCK_HEADER_LENGTH = 18;
const int BLOCK_FOOTER_LENGTH = 8;

// compresses the current block
int CompressBlock(char* uncompressedBlock, unsigned int* pUncompressLen, char* compressedBlock, const unsigned int compressedLen, const int compressionLevel);
// augmented the gzread parameters to use offsets rather than pointers
int gzreadOffset(gzFile file, char* buf, int offset, unsigned len);
// augmented the gzwrite parameters to use offsets rather than pointers
int gzwriteOffset(gzFile file, char* buf, int offset, unsigned len);
// packs an unsigned integer into the specified buffer
void PackUnsignedInt(char* buffer, const unsigned int value);
// packs an unsigned short into the specified buffer
void PackUnsignedShort(char* buffer, unsigned short value);
// uncompresses the current block
int UncompressBlock(const char* compressedBlock, const unsigned int compressedLen, char* uncompressedBlock, const unsigned int uncompressedLen);
