/*
 * tga_writer.c
 *
 *  Created on: Apr 27, 2013
 *      Author: Sachith Dunatunga
 */
#include <stdio.h>

#include "tga_writer.h"

//#ifdef __cplusplus
//extern "C" {
//#endif

void create_tga_header(tgaheader_t *th, int width, int height, int bitdepth)
{
	/* Uncompressed header:
	 * {0,0,2,0,0,0,0,0,XPIXELS,YPIXELS,24,0,NULL}; */
	th->idlength = 0;
	th->colourmaptype = 0;
	th->datatypecode = 2;
	th->colourmaporigin = 0;
	th->colourmaplength = 0;
	th->colourmapdepth = 0;
	th->x_origin = 0;
	th->y_origin = 0;
	th->width = width;
	th->height = height;
	th->bitsperpixel = bitdepth;
	th->imagedescriptor = 0;
	th->idfield[0] = 0;

	return;
}

void write_tga_header(FILE *fp, tgaheader_t *th)
{
	int i;

	putc(th->idlength, fp);
	putc(th->colourmaptype, fp);
	putc(th->datatypecode, fp);

	putc((th->colourmaporigin & 0x00FF), fp);
	putc((th->colourmaporigin & 0xFF00) >> 8, fp);

	putc((th->colourmaplength & 0x00FF), fp);
	putc((th->colourmaplength & 0xFF00) >> 8, fp);

	putc(th->colourmapdepth, fp);

	putc((th->x_origin & 0x00FF), fp);
	putc((th->x_origin & 0xFF00) >> 8, fp);
	putc((th->y_origin & 0x00FF), fp);
	putc((th->y_origin & 0xFF00) >> 8, fp);

	putc((th->width & 0x00FF), fp);
	putc((th->width & 0xFF00) >> 8, fp);
	putc((th->height & 0x00FF), fp);
	putc((th->height & 0xFF00) >> 8, fp);

	putc(th->bitsperpixel, fp);
	putc(th->imagedescriptor, fp);

	for (i = 0; i < th->idlength; i++) {
		putc(th->idfield[i], fp);
	}
	return;
}

void write_tga_image(FILE *fp, tgaheader_t *th, char *tga_pixeldata)
{
	int i;
	int pixeldata_len;

	pixeldata_len = (th->width * th->height * th->bitsperpixel) >> 3;
	write_tga_header(fp, th);
	for (i = 0; i < pixeldata_len; i++) {
		putc(tga_pixeldata[i], fp);
	}
	return;
}

//#ifdef __cplusplus
//}
//#endif
