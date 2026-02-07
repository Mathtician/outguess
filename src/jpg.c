/*
 * Copyright 1999-2001 Niels Provos <provos@citi.umich.edu>
 * Copyright 2016      Joao Eriberto Mota Filho <eriberto@debian.org>
 * Copyright 2020      Robin Vobruba <hoijui.quaero@gmail.com>
 * Copyright 2021      Daniel T. Borelli <daltomi@disroot.org>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *      This product includes software developed by Niels Provos.
 * 4. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * This file is derived from example.c distributed with the jpeg-6b
 * source distribution. See the JPEG-README file.
 */

#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "config.h"
#include "outguess.h"
#include "pnm.h"
#include <jpeglib.h>
#include <jmorecfg.h>
#include "jpg.h"

struct og_jpeg_marker {
	int marker;
	unsigned int data_len;
	unsigned char *data;
	struct og_jpeg_marker *next;
};

static int
jpeg_marker_is_icc(jpeg_saved_marker_ptr marker)
{
	static const unsigned char icc_id[] = "ICC_PROFILE";
	const size_t icc_id_len = sizeof(icc_id) - 1;

	if (marker->marker != (JPEG_APP0 + 2))
		return 0;
	if (marker->data_length < icc_id_len + 3)
		return 0;
	if (memcmp(marker->data, icc_id, icc_id_len) != 0)
		return 0;
	return marker->data[icc_id_len] == '\0';
}

static struct og_jpeg_marker *
jpeg_markers_copy_icc(struct jpeg_decompress_struct *cinfo)
{
	struct og_jpeg_marker *head = NULL;
	struct og_jpeg_marker **tail = &head;
	jpeg_saved_marker_ptr marker;

	for (marker = cinfo->marker_list; marker != NULL;
	        marker = marker->next) {
		struct og_jpeg_marker *node;

		if (!jpeg_marker_is_icc(marker))
			continue;

		node = checkedmalloc(sizeof(*node));
		node->marker = marker->marker;
		node->data_len = marker->data_length;
		node->data = checkedmalloc(node->data_len);
		memcpy(node->data, marker->data, node->data_len);
		node->next = NULL;

		*tail = node;
		tail = &node->next;
	}

	return head;
}

static void
jpeg_markers_write(
        struct jpeg_compress_struct *cinfo, struct og_jpeg_marker *markers)
{
	struct og_jpeg_marker *marker;

	for (marker = markers; marker != NULL; marker = marker->next)
		jpeg_write_marker(
		        cinfo, marker->marker, marker->data, marker->data_len);
}

void
jpeg_markers_free(struct og_jpeg_marker *markers)
{
	while (markers != NULL) {
		struct og_jpeg_marker *next = markers->next;

		free(markers->data);
		free(markers);
		markers = next;
	}
}

void jpeg_dummy_dest(j_compress_ptr cinfo);

static void jpeg_compress_image(FILE *outfile, image *image);
static bitmap *jpeg_bitmap_from_coefficients(FILE *infile, int eval);
static void jpeg_traverse_coeffs(
        j_decompress_ptr cinfo, jvirt_barray_ptr *coef_arrays, int write_mode);

/* The functions that can be used to handle a JPEG data object */

handler jpg_handler = {"jpg", NULL, init_JPEG_handler, read_JPEG_file,
        write_JPEG_file, bitmap_from_jpg, bitmap_to_jpg, preserve_jpg};

static int jpeg_state;
static bitmap tbitmap;
static u_int32_t off;
static int quality = 75;
static int jpeg_eval;
static int eval_cnt;

static int dctmin;
static int dctmax;

extern int steg_foil; /* Statistics keps in main program */
extern int steg_foilfail;

#define DCTMIN 100
#define DCTENTRIES 256
static int dctadjust[DCTENTRIES];

#define DCTFREQRANGE 5000 /* Number of bits for which the below holds */
#define DCTFREQREDUCE 33  /* Threshold is /REDUCE, 1% for 100 */
#define DCTFREQMIN 2      /* At least 5 coeff in cache */
static int dctfreq[DCTENTRIES];
static int dctpending;

void
init_state(int state, int eval, bitmap *bitmap)
{
	jpeg_state = state;
	jpeg_eval = eval;
	eval_cnt = 0;

	dctmin = 127;
	dctmax = -127;

	off = 0;
	if (state == JPEG_READING) {
		memset(&tbitmap, 0, sizeof(tbitmap));
		tbitmap.bytes = 256;
		tbitmap.bits = tbitmap.bytes * 8;
		tbitmap.bitmap = checkedmalloc(tbitmap.bytes);
		tbitmap.locked = checkedmalloc(tbitmap.bytes);
		memset(tbitmap.locked, 0, tbitmap.bytes);
		tbitmap.data = checkedmalloc(tbitmap.bits);
	} else if (bitmap) {
		memcpy(&tbitmap, bitmap, sizeof(tbitmap));
	}
}

int
preserve_single(bitmap *bitmap, int off, char coeff)
{
	int i;
	char *data = bitmap->data;
	char *pbits = bitmap->bitmap;
	char *plock = bitmap->locked;
	char *pmetalock = bitmap->metalock;

	for (i = off - 1; i >= 0; i--) {
		if (TEST_BIT(plock, i))
			continue;
		if (TEST_BIT(pmetalock, i))
			continue;

		/* Switch the coefficient to the value that we just replaced */
		if (data[i] == coeff) {
			char cbit;

			data[i] = coeff ^ 0x01;

			cbit = (unsigned char)coeff & 0x01;
			WRITE_BIT(pbits, i, cbit ^ 0x01);

			WRITE_BIT(pmetalock, i, 1);

			if (jpeg_eval)
				fprintf(stderr,
				        "off: %d, i: %d, coeff: %d, data: %d\n",
				        off, i, coeff, data[i]);

			return (i);
		}
	}

	return (-1);
}

int
preserve_jpg(bitmap *bitmap, int off)
{
	char coeff;
	char *data = bitmap->data;

	if (off == -1) {
		int res;

		if (jpeg_eval)
			fprintf(stderr, "DCT: %d<->%d\n", dctmin, dctmax);

		bitmap->preserve = preserve_jpg;
		memset(bitmap->metalock, 0, bitmap->bytes);

		memset(dctadjust, 0, sizeof(dctadjust));
		memset(dctfreq, 0, sizeof(dctfreq));
		dctpending = 0;

		/* Calculate coefficent frequencies */
		for (int i = 0; i < bitmap->bits; i++) {
			dctfreq[data[i] + 127]++;
		}

		int a = dctfreq[-1 + 127];
		int b = dctfreq[-2 + 127];

		if (a < b) {
			fprintf(stderr, "Can not calculate estimate\n");
			res = -1;
		} else
			res = 2 * bitmap->bits * b / (a + b);

		/* Pending threshold based on frequencies */
		for (int i = 0; i < DCTENTRIES; i++) {
			dctfreq[i] = dctfreq[i] /
			             ((float)bitmap->bits / DCTFREQRANGE);
			dctfreq[i] /= DCTFREQREDUCE;
			if (dctfreq[i] < DCTFREQMIN)
				dctfreq[i] = DCTFREQMIN;

			if (jpeg_eval)
				fprintf(stderr, "Foil: %d :< %d\n", i - 127,
				        dctfreq[i]);
		}

		bitmap->maxcorrect = (res < 0) ? 0 : (size_t)res;
		return (res);
	} else if (off >= bitmap->bits) {
		/* Reached end of image */
		for (int i = 0; i < DCTENTRIES; i++) {
			while (dctadjust[i]) {
				dctadjust[i]--;

				coeff = i - 127;

				if (preserve_single(bitmap, bitmap->bits - 1,
				            coeff) != -1)
					steg_foil++;
				else
					steg_foilfail++;
			}
		}

		return (0);
	}

	/* We need to find this coefficient, and change it to data[off] */
	coeff = data[off] ^ 0x01;

	if (dctadjust[data[off] + 127]) {
		/* But we are still missing compensation for the opposite */
		dctadjust[data[off] + 127]--;
		dctpending--;
		return (0);
	}

	if (dctadjust[coeff + 127] < dctfreq[coeff + 127]) {
		dctadjust[coeff + 127]++;
		dctpending++;
		return (0);
	}

	const int i = preserve_single(bitmap, off, coeff);

	if (i != -1) {
		steg_foil++;
		return (i);
	}

	/* We have one too many of this */
	dctadjust[coeff + 127]++;
	dctpending++;

	return (-1);
}

bitmap *
finish_state(void)
{
	int i;
	bitmap *pbitmap;

	if (jpeg_eval)
		fprintf(stderr, "\n");

	if (jpeg_state != JPEG_READING)
		return NULL;

	tbitmap.bits = off;
	tbitmap.bytes = (off + 7) / 8;

	tbitmap.detect = checkedmalloc(tbitmap.bits);
	tbitmap.metalock = checkedmalloc(tbitmap.bytes);

	for (i = 0; i < off; i++) {
		char temp = abs(tbitmap.data[i]);
		if (temp >= JPG_THRES_MAX)
			tbitmap.detect[i] = -1;
		else if (temp >= JPG_THRES_LOW)
			tbitmap.detect[i] = 0;
		else if (temp >= JPG_THRES_MIN)
			tbitmap.detect[i] = 1;
		else
			tbitmap.detect[i] = 2;
	}

	pbitmap = checkedmalloc(sizeof(bitmap));

	memcpy(pbitmap, &tbitmap, sizeof(tbitmap));

	return pbitmap;
}

short
steg_use_bit(unsigned short temp)
{
	if (jpeg_state == JPEG_DISABLED)
		return temp;

	if ((temp & 0x1) == temp)
		goto steg_end;

	switch (jpeg_state) {
	case JPEG_READING:
		WRITE_BIT(tbitmap.bitmap, off, temp & 0x1);
		tbitmap.data[off] = temp;

		if ((short)temp < dctmin)
			dctmin = (short)temp;
		if ((short)temp > dctmax)
			dctmax = (short)temp;

		off++;

		if (off >= tbitmap.bits) {
			u_char *buf;

			tbitmap.bytes += 256;
			tbitmap.bits += 256 * 8;
			if (!(buf = realloc(tbitmap.bitmap, tbitmap.bytes))) {
				fprintf(stderr, "steg_use_bit: realloc()\n");
				exit(1);
			}
			tbitmap.bitmap = buf;
			if (!(buf = realloc(tbitmap.locked, tbitmap.bytes))) {
				fprintf(stderr, "steg_use_bit: realloc()\n");
				exit(1);
			}
			tbitmap.locked = buf;
			memset(tbitmap.locked + tbitmap.bytes - 256, 0, 256);
			if (!(buf = realloc(tbitmap.data, tbitmap.bits))) {
				fprintf(stderr, "steg_use_bit: realloc()\n");
				exit(1);
			}
			tbitmap.data = buf;
		}
		break;
	default:
		temp = (temp & ~0x1) | (TEST_BIT(tbitmap.bitmap, off) ? 1 : 0);
		off++;

		break;
	}

steg_end:
	if (jpeg_eval) {
		if (eval_cnt % DCTSIZE2 == 0)
			fprintf(stderr, "\n[%d]%.7d: ", jpeg_state, eval_cnt);
		if ((temp & 0x1) != temp)
			fprintf(stderr, "% .3d,", (short)temp);
		eval_cnt++;
	}

	return temp;
}

void
init_JPEG_handler(char *parameter)
{
	if (parameter)
		quality = atoi(parameter);
	if (quality < 75)
		quality = 75;
	if (quality > 100)
		quality = 100;
	fprintf(stderr, "JPEG compression quality set to %d\n", quality);
}

void
bitmap_from_jpg(bitmap *dbitmap, image *image, int flags)
{
	bitmap *tmpmap;

	if (flags & STEG_RETRIEVE) {
		memcpy(dbitmap, image->bitmap, sizeof(*dbitmap));
		free(image->bitmap);
		image->bitmap = NULL;
		return;
	}

	if (image->bitmap) {
		tmpmap = image->bitmap;
		free(tmpmap->bitmap);
		free(tmpmap->locked);
		free(tmpmap);
		image->bitmap = NULL;
	}

	tmpmap = compress_JPEG(image);
	memcpy(dbitmap, tmpmap, sizeof(*dbitmap));
	free(tmpmap);
}

void
bitmap_to_jpg(image *image, bitmap *bitmap, int flags)
{
	init_state(JPEG_WRITING, steg_stat >= 3 ? 1 : 0, bitmap);
}

/******************** JPEG COMPRESSION SAMPLE INTERFACE *******************/

/* This half of the example shows how to feed data into the JPEG compressor.
 * We present a minimal version that does not worry about refinements such
 * as error recovery (the JPEG code will just exit() if it gets an error).
 */

/*
 * IMAGE DATA FORMATS:
 *
 * The standard input image format is a rectangular array of pixels, with
 * each pixel having the same number of "component" values (color channels).
 * Each pixel row is an array of JSAMPLEs (which typically are unsigned chars).
 * If you are working with color data, then the color values for each pixel
 * must be adjacent in the row; for example, R,G,B,R,G,B,R,G,B,... for 24-bit
 * RGB color.
 *
 * For this example, we'll assume that this data structure matches the way
 * our application has stored the image in memory, so we can just pass a
 * pointer to our image buffer.  In particular, let's say that the image is
 * RGB color and is described by:
 */

extern JSAMPLE *image_buffer; /* Points to large array of R,G,B-order data */
extern int image_height;      /* Number of rows in image */
extern int image_width;       /* Number of columns in image */

bitmap *
compress_JPEG(image *image)
{
	FILE *tmp;
	bitmap *bmp;

	tmp = tmpfile();
	if (tmp == NULL) {
		fprintf(stderr, "compress_JPEG: tmpfile()\n");
		exit(1);
	}

	jpeg_compress_image(tmp, image);
	fflush(tmp);
	rewind(tmp);

	bmp = jpeg_bitmap_from_coefficients(tmp, steg_stat >= 3 ? 1 : 0);
	fclose(tmp);

	return bmp;
}

/*
 * Sample routine for JPEG compression.  We assume that the target file name
 * and a compression quality factor are passed in.
 */

void
write_JPEG_file(FILE *outfile, image *image)
{
	struct jpeg_decompress_struct srcinfo;
	struct jpeg_compress_struct dstinfo;
	struct jpeg_error_mgr jsrcerr;
	struct jpeg_error_mgr jdsterr;
	jvirt_barray_ptr *coef_arrays;
	FILE *tmp;

	tmp = tmpfile();
	if (tmp == NULL) {
		fprintf(stderr, "write_JPEG_file: tmpfile()\n");
		exit(1);
	}

	jpeg_compress_image(tmp, image);
	fflush(tmp);
	rewind(tmp);

	srcinfo.err = jpeg_std_error(&jsrcerr);
	jpeg_create_decompress(&srcinfo);
	jpeg_stdio_src(&srcinfo, tmp);
	(void)jpeg_read_header(&srcinfo, TRUE);

	coef_arrays = jpeg_read_coefficients(&srcinfo);
	if (coef_arrays == NULL) {
		fprintf(stderr, "write_JPEG_file: jpeg_read_coefficients()\n");
		exit(1);
	}

	/* Apply the embedding bits directly to the quantized DCT coefficients. */
	jpeg_traverse_coeffs(&srcinfo, coef_arrays, 1);

	dstinfo.err = jpeg_std_error(&jdsterr);
	jpeg_create_compress(&dstinfo);
	jpeg_stdio_dest(&dstinfo, outfile);
	jpeg_copy_critical_parameters(&srcinfo, &dstinfo);
	jpeg_write_coefficients(&dstinfo, coef_arrays);
	if (image->jpeg_markers != NULL)
		jpeg_markers_write(&dstinfo, image->jpeg_markers);

	jpeg_finish_compress(&dstinfo);
	jpeg_destroy_compress(&dstinfo);

	(void)jpeg_finish_decompress(&srcinfo);
	jpeg_destroy_decompress(&srcinfo);
	fclose(tmp);

	/* After finish_compress, we can close the output file. */
	fclose(outfile);
}

/*
 * SOME FINE POINTS:
 *
 * In the above loop, we ignored the return value of jpeg_write_scanlines,
 * which is the number of scanlines actually written.  We could get away
 * with this because we were only relying on the value of cinfo.next_scanline,
 * which will be incremented correctly.  If you maintain additional loop
 * variables then you should be careful to increment them properly.
 * Actually, for output to a stdio stream you needn't worry, because
 * then jpeg_write_scanlines will write all the lines passed (or else exit
 * with a fatal error).  Partial writes can only occur if you use a data
 * destination module that can demand suspension of the compressor.
 * (If you don't know what that's for, you don't need it.)
 *
 * If the compressor requires full-image buffers (for entropy-coding
 * optimization or a multi-scan JPEG file), it will create temporary
 * files for anything that doesn't fit within the maximum-memory setting.
 * (Note that temp files are NOT needed if you use the default parameters.)
 * On some systems you may need to set up a signal handler to ensure that
 * temporary files are deleted if the program is interrupted.  See libjpeg.doc.
 *
 * Scanlines MUST be supplied in top-to-bottom order if you want your JPEG
 * files to be compatible with everyone else's.  If you cannot readily read
 * your data in that order, you'll need an intermediate array to hold the
 * image.  See rdtarga.c or rdbmp.c for examples of handling bottom-to-top
 * source data using the JPEG code's internal virtual-array mechanisms.
 */

/******************** JPEG DECOMPRESSION SAMPLE INTERFACE *******************/

/* This half of the example shows how to read data from the JPEG decompressor.
 * It's a bit more refined than the above, in that we show:
 *   (a) how to modify the JPEG library's standard error-reporting behavior;
 *   (b) how to allocate workspace using the library's memory manager.
 *
 * Just to make this example a little different from the first one, we'll
 * assume that we do not intend to put the whole image into an in-memory
 * buffer, but to send it line-by-line someplace else.  We need a one-
 * scanline-high JSAMPLE array as a work buffer, and we will let the JPEG
 * memory manager allocate it for us.  This approach is actually quite useful
 * because we don't need to remember to deallocate the buffer separately: it
 * will go away automatically when the JPEG object is cleaned up.
 */

/*
 * Sample routine for JPEG decompression.  We assume that the source file name
 * is passed in.  We want to return 1 on success, 0 on error.
 */

image *
read_JPEG_file(FILE *infile)
{
	/* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
	struct jpeg_decompress_struct cinfo;
	/* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
	struct jpeg_error_mgr jerr;

	/* More stuff */
	image *image;
	JSAMPARRAY buffer; /* Output row buffer */
	int row_stride;    /* physical row width in output buffer */
	int saved_state;
	FILE *src;
	FILE *tmp;
	int seekable;

	image = checkedmalloc(sizeof(*image));
	memset(image, 0, sizeof(*image));

	src = infile;
	tmp = NULL;
	seekable = (fseek(infile, 0, SEEK_CUR) == 0);
	if (!seekable) {
		char buf[8192];
		size_t nread;

		tmp = tmpfile();
		if (tmp == NULL) {
			fprintf(stderr, "read_JPEG_file: tmpfile()\n");
			exit(1);
		}

		while ((nread = fread(buf, 1, sizeof(buf), infile)) > 0) {
			if (fwrite(buf, 1, nread, tmp) != nread) {
				fprintf(stderr, "read_JPEG_file: fwrite()\n");
				exit(1);
			}
		}
		if (ferror(infile)) {
			fprintf(stderr, "read_JPEG_file: fread()\n");
			exit(1);
		}

		src = tmp;
	}

	if (fseek(src, 0, SEEK_SET) != 0) {
		fprintf(stderr, "read_JPEG_file: fseek()\n");
		exit(1);
	}

	image->bitmap = jpeg_bitmap_from_coefficients(src, 0);
	if (fseek(src, 0, SEEK_SET) != 0) {
		fprintf(stderr, "read_JPEG_file: fseek()\n");
		exit(1);
	}

	/* Step 1: allocate and initialize JPEG decompression object */

	cinfo.err = jpeg_std_error(&jerr);
	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_decompress(&cinfo);

	/* Step 2: specify data source (eg, a file) */

	jpeg_stdio_src(&cinfo, src);

	/* Step 3: read file parameters with jpeg_read_header() */

	jpeg_save_markers(&cinfo, JPEG_APP0 + 2, 0xFFFF);
	(void)jpeg_read_header(&cinfo, TRUE);
	image->jpeg_markers = jpeg_markers_copy_icc(&cinfo);
	/* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.doc for more info.
   */

	/* Step 4: set parameters for decompression */

	/* In this example, we don't need to change any of the defaults set by
   * jpeg_read_header(), so we do nothing here.
   */

	/* Step 5: Start decompressor */

	saved_state = jpeg_state;
	jpeg_state = JPEG_DISABLED;

	(void)jpeg_start_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

	/* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */

	image->x = cinfo.output_width;
	image->y = cinfo.output_height;
	image->depth = cinfo.output_components;
	image->max = 255;

	image->img = checkedmalloc(cinfo.output_width * cinfo.output_height *
	                           cinfo.output_components);

	/* JSAMPLEs per row in output buffer */
	row_stride = cinfo.output_width * cinfo.output_components;
	/* Make a one-row-high sample array that will go away when done with image */
	buffer = (*cinfo.mem->alloc_sarray)(
	        (j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);

	/* Step 6: while (scan lines remain to be read) */
	/*           jpeg_read_scanlines(...); */

	/* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */
	while (cinfo.output_scanline < cinfo.output_height) {
		/* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
		(void)jpeg_read_scanlines(&cinfo, buffer, 1);
		/* Assume put_scanline_someplace wants a pointer and sample count. */

		memcpy(&image->img[(cinfo.output_scanline - 1) * row_stride],
		        buffer[0], row_stride);
	}

	/* Step 7: Finish decompression */

	(void)jpeg_finish_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */
	jpeg_state = saved_state;

	/* Step 8: Release JPEG decompression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_decompress(&cinfo);

	if (tmp != NULL)
		fclose(tmp);

	/* And we're done! */
	return image;
}

/*
 * SOME FINE POINTS:
 *
 * In the above code, we ignored the return value of jpeg_read_scanlines,
 * which is the number of scanlines actually read.  We could get away with
 * this because we asked for only one line at a time and we weren't using
 * a suspending data source.  See libjpeg.doc for more info.
 *
 * We cheated a bit by calling alloc_sarray() after jpeg_start_decompress();
 * we should have done it beforehand to ensure that the space would be
 * counted against the JPEG max_memory setting.  In some systems the above
 * code would risk an out-of-memory error.  However, in general we don't
 * know the output image dimensions before jpeg_start_decompress(), unless we
 * call jpeg_calc_output_dimensions().  See libjpeg.doc for more about this.
 *
 * Scanlines are returned in the same order as they appear in the JPEG file,
 * which is standardly top-to-bottom.  If you must emit data bottom-to-top,
 * you can use one of the virtual arrays provided by the JPEG memory manager
 * to invert the data.  See wrbmp.c for an example.
 *
 * As with compression, some operating modes may require temporary files.
 * On some systems you may need to set up a signal handler to ensure that
 * temporary files are deleted if the program is interrupted.  See libjpeg.doc.
 */

/* Expanded data destination object for dummy output */

typedef struct {
	struct jpeg_destination_mgr pub; /* public fields */
} my_destination_mgr;

#define BUFSIZE 256

char dummy_buf[BUFSIZE];

typedef my_destination_mgr *my_dest_ptr;

METHODDEF(void)
init_destination(j_compress_ptr cinfo)
{
	my_dest_ptr dest = (my_dest_ptr)cinfo->dest;

	dest->pub.next_output_byte = dummy_buf;
	dest->pub.free_in_buffer = BUFSIZE;
}

METHODDEF(boolean)
empty_output_buffer(j_compress_ptr cinfo)
{
	my_dest_ptr dest = (my_dest_ptr)cinfo->dest;

	dest->pub.next_output_byte = dummy_buf;
	dest->pub.free_in_buffer = BUFSIZE;

	return TRUE;
}

METHODDEF(void)
term_destination(j_compress_ptr cinfo)
{
}

void
jpeg_dummy_dest(j_compress_ptr cinfo)
{
	my_dest_ptr dest;

	/* The destination object is made permanent so that multiple JPEG images
   * can be written to the same file without re-executing jpeg_stdio_dest.
   * This makes it dangerous to use this manager and a different destination
   * manager serially with the same JPEG object, because their private object
   * sizes may be different.  Caveat programmer.
   */
	if (cinfo->dest == NULL) { /* first time for this JPEG object? */
		cinfo->dest = (struct jpeg_destination_mgr *)(*cinfo->mem
		                ->alloc_small)((j_common_ptr)cinfo,
		        JPOOL_PERMANENT, sizeof(my_destination_mgr));
	}

	dest = (my_dest_ptr)cinfo->dest;
	dest->pub.init_destination = init_destination;
	dest->pub.empty_output_buffer = empty_output_buffer;
	dest->pub.term_destination = term_destination;
}

static void
jpeg_traverse_coeffs(
        j_decompress_ptr cinfo, jvirt_barray_ptr *coef_arrays, int write_mode)
{
	JDIMENSION mcu_row, mcu_col;

	for (mcu_row = 0; mcu_row < cinfo->total_iMCU_rows; mcu_row++) {
		for (mcu_col = 0; mcu_col < cinfo->MCUs_per_row; mcu_col++) {
			int ci;

			for (ci = 0; ci < cinfo->num_components; ci++) {
				jpeg_component_info *compptr =
				        cinfo->comp_info + ci;
				int mcu_width = compptr->MCU_width;
				int mcu_height = compptr->MCU_height;
				int useful_width =
				        (mcu_col < cinfo->MCUs_per_row - 1)
				                ? mcu_width
				                : compptr->last_col_width;
				int useful_height =
				        (mcu_row < cinfo->total_iMCU_rows - 1)
				                ? mcu_height
				                : compptr->last_row_height;
				int y, x, k;

				for (y = 0; y < useful_height; y++) {
					JDIMENSION block_row =
					        mcu_row * mcu_height + y;
					JBLOCKARRAY buffer = (*cinfo->mem
					                ->access_virt_barray)(
					        (j_common_ptr)cinfo,
					        coef_arrays[ci], block_row, 1,
					        write_mode ? TRUE : FALSE);
					JBLOCKROW row = buffer[0];
					JDIMENSION block_col_base =
					        mcu_col * mcu_width;

					for (x = 0; x < useful_width; x++) {
						JCOEFPTR block =
						        row[block_col_base + x];

						for (k = 0; k < DCTSIZE2; k++) {
							if (write_mode)
								block[k] = (JCOEF)steg_use_bit(
								        (unsigned short)block
								                [k]);
							else
								(void)steg_use_bit(
								        (unsigned short)block
								                [k]);
						}
					}
				}
			}
		}
	}
}

static bitmap *
jpeg_bitmap_from_coefficients(FILE *infile, int eval)
{
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	jvirt_barray_ptr *coef_arrays;
	bitmap *bmp;

	init_state(JPEG_READING, eval, NULL);

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, infile);
	(void)jpeg_read_header(&cinfo, TRUE);

	coef_arrays = jpeg_read_coefficients(&cinfo);
	if (coef_arrays == NULL) {
		fprintf(stderr, "jpeg_bitmap_from_coefficients: "
		                "jpeg_read_coefficients()\n");
		exit(1);
	}

	jpeg_traverse_coeffs(&cinfo, coef_arrays, 0);
	bmp = finish_state();

	(void)jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);

	return bmp;
}

static void
jpeg_compress_image(FILE *outfile, image *image)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW row_pointer[1];
	int row_stride;
	int saved_state;

	saved_state = jpeg_state;
	jpeg_state = JPEG_DISABLED;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, outfile);

	cinfo.image_width = image->x;
	cinfo.image_height = image->y;
	cinfo.input_components = image->depth;
	cinfo.in_color_space = JCS_RGB;

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(
	        &cinfo, quality, TRUE /* limit to baseline-JPEG values */);
	jpeg_start_compress(&cinfo, TRUE);

	row_stride = image->x * 3;
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &image->img[cinfo.next_scanline * row_stride];
		(void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	jpeg_state = saved_state;
}
