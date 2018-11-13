/* Wrappers around zlib to make interfacing to it a bit easier. */

/* Copyright (C) 2009 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "./libdeflate/libdeflate.h"

#define DEFLATE_COMPRESSION_LEVEL 6

/**
 * Convert error code to message
 */
static const char *libdeflate_result_to_string(enum libdeflate_result result) {
    switch (result) {
        /* Decompression was successful.  */
        case LIBDEFLATE_SUCCESS: return "success";

        /* Decompressed failed because the compressed data was invalid, corrupt,
         * or otherwise unsupported.  */
        case LIBDEFLATE_BAD_DATA: return "corrupt or invalid data";

        /* A NULL 'actual_out_nbytes_ret' was provided, but the data would have
         * decompressed to fewer than 'out_nbytes_avail' bytes.  */
        case LIBDEFLATE_SHORT_OUTPUT: return "data return size is shorter than expected";

        /* The data would have decompressed to more than 'out_nbytes_avail'
         * bytes.  */
        case LIBDEFLATE_INSUFFICIENT_SPACE: return "not enough space to decompress data";
    }

    return "unknown result";
}

/**
 * Return size of buffer needed to compress something of given size uncompressed.
 */
size_t zCompBufSize(size_t uncompressedSize)
{
    return 1.001 * uncompressedSize + 13;
}

/**
 * Compress data from memory to memory.  Returns size after compression.
 */
size_t zCompress(
        void *uncompressed,       /* Start of area to compress. */
        size_t uncompressedSize,  /* Size of area to compress. */
        void *compressed,         /* Where to put compressed bits */
        size_t compressedSize)    /* Size of compressed bits - calculate using zCompBufSize */
{
    static struct libdeflate_compressor *compressor = NULL;

    if (compressor == NULL)
        compressor = libdeflate_alloc_compressor(DEFLATE_COMPRESSION_LEVEL);

    size_t size = libdeflate_zlib_compress(compressor, uncompressed,
                                           uncompressedSize, compressed, compressedSize);

    if (size == 0)
        errAbort("Error: couldn't compress data");

    return size;
}

/*
 * Uncompress data from memory to memory.  Returns size after decompression.
 */
size_t zUncompress(
        void *compressed,       /* Compressed area */
        size_t compressedSize,  /* Size after compression */
        void *uncompressed,        /* Where to put uncompressed bits */
        size_t uncompressedSize)   /* Max size of uncompressed bits. */
{
    static struct libdeflate_decompressor *decompressor = NULL;

    if (decompressor == NULL)
        decompressor = libdeflate_alloc_decompressor();

    size_t actualSize = 0;
    enum libdeflate_result result = libdeflate_zlib_decompress(decompressor, compressed, compressedSize,
                                                               uncompressed, uncompressedSize, &actualSize);

    if (result != 0)
        errAbort("Error: couldn't decompress data: %s", libdeflate_result_to_string(result));

    return actualSize;
}

void zSelfTest(int count)
/* Run an internal diagnostic. */
{
bits32 testData[count];
int uncSize = count*sizeof(bits32);
int i;
for (i=0; i<count; ++i)
    testData[i] = i;
int compressedSize = zCompBufSize(uncSize);
char compressed[compressedSize];
int compSize = zCompress(testData, uncSize, compressed, compressedSize);
char uncBuf[uncSize];
zUncompress(compressed, compSize, uncBuf, uncSize);
if (memcmp(uncBuf, testData, uncSize) != 0)
    errAbort("zSelfTest %d failed", count);
else
    verbose(2, "zSelfTest %d passed, compression ratio %3.1f\n", count, (double)compSize/uncSize);
}
