#pragma once

#include "BGZF.h"

// compresses the current block
int CompressBlock(char* uncompressedBlock, unsigned int* pUncompressLen, char* compressedBlock, const unsigned int compressedLen, const int compressionLevel) {

	// initialization
	char* buffer = compressedBlock;

	int inputLength      = *pUncompressLen;
	int compressedLength = 0;
	int status;
	int remaining;

	unsigned int bufferSize = compressedLen;
	unsigned int crc;

	z_stream zs;

	// initialize the gzip header
	memset(buffer, 0, 18);
	buffer[0]  = GZIP_ID1;
	buffer[1]  = (char)GZIP_ID2;
	buffer[2]  = CM_DEFLATE;
	buffer[3]  = FLG_FEXTRA;
	buffer[9]  = (char)OS_UNKNOWN;
	buffer[10] = BGZF_XLEN;
	buffer[12] = BGZF_ID1;
	buffer[13] = BGZF_ID2;
	buffer[14] = BGZF_LEN;

	// loop to retry for blocks that do not compress enough
	while(1) {

		// initialize zstream values
		zs.zalloc    = NULL;
		zs.zfree     = NULL;
		zs.next_in   = (Bytef*)uncompressedBlock;
		zs.avail_in  = inputLength;
		zs.next_out  = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
		zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

		// initialize the zlib compression algorithm
		if(deflateInit2(&zs, compressionLevel, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK) {
			printf("ERROR: zlib deflate initialization failed.\n");
			return -1;
		}

		// compress the data
		status = deflate(&zs, Z_FINISH);
		if(status != Z_STREAM_END) {

			deflateEnd(&zs);

			// reduce the input length and try again
			if(status == Z_OK) {
				inputLength -= 1024;
				if(inputLength < 0) {
					printf("ERROR: input reduction failed.\n");
					return -1;
				}
				continue;
			}

			printf("ERROR: zlib deflate failed.\n");
			return -1;
		}

		// finalize the compression routine
		if(deflateEnd(&zs) != Z_OK) {
			printf("ERROR: deflate end failed.\n");
			return -1;
		}

		compressedLength = zs.total_out;
		compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

		if(compressedLength > MAX_BLOCK_SIZE) {
			printf("ERROR: deflate overflow.\n");
			return -1;
		}

		break;
	}

	// store the compressed length
	PackUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

	// store the CRC32 checksum
	crc = crc32(0, NULL, 0);
	crc = crc32(crc, (Bytef*)uncompressedBlock, inputLength);
	PackUnsignedInt(&buffer[compressedLength - 8], crc);
	PackUnsignedInt(&buffer[compressedLength - 4], inputLength);

	// ensure that we have less than a block of data left
	remaining = *pUncompressLen - inputLength;
	if(remaining > 0) {
		if(remaining > inputLength) {
			printf("ERROR: remainder too large.\n");
			return -1;
		}
		memcpy(uncompressedBlock, uncompressedBlock + inputLength, remaining);
	}

	*pUncompressLen = remaining;
	return compressedLength;
}

// augmented the gzread parameters to use offsets rather than pointers
int gzreadOffset(gzFile file, char* buf, int offset, unsigned len) {
	return gzread(file, buf + offset, len);
}

// augmented the gzwrite parameters to use offsets rather than pointers
int gzwriteOffset(gzFile file, char* buf, int offset, unsigned len) {
	return gzwrite(file, buf + offset, len);
}

// packs an unsigned integer into the specified buffer
void PackUnsignedInt(char* buffer, const unsigned int value) {
	buffer[0] = (char)value;
	buffer[1] = (char)(value >> 8);
	buffer[2] = (char)(value >> 16);
	buffer[3] = (char)(value >> 24);
}

// packs an unsigned short into the specified buffer
void PackUnsignedShort(char* buffer, unsigned short value) {
	buffer[0] = (char)value;
	buffer[1] = (char)(value >> 8);
}

// uncompresses the current block
int UncompressBlock(const char* compressedBlock, const unsigned int compressedLen, char* uncompressedBlock, const unsigned int uncompressedLen) {

	// initialization
	z_stream zs;
	int status;

	// Inflate the block in m_BGZF.CompressedBlock into m_BGZF.UncompressedBlock
	zs.zalloc    = NULL;
	zs.zfree     = NULL;
	zs.next_in   = (Bytef*)compressedBlock + 18;
	zs.avail_in  = compressedLen - 16;
	zs.next_out  = (Bytef*)uncompressedBlock;
	zs.avail_out = uncompressedLen;

	status = inflateInit2(&zs, GZIP_WINDOW_BITS);
	if (status != Z_OK) {
		printf("inflateInit failed\n");
		return -1;
	}

	status = inflate(&zs, Z_FINISH);
	if (status != Z_STREAM_END) {
		inflateEnd(&zs);
		printf("ERROR: inflate failed\n");
		return -1;
	}

	status = inflateEnd(&zs);
	if (status != Z_OK) {
		printf("ERROR: nflateEnd failed\n");
		return -1;
	}

	return zs.total_out;
}
