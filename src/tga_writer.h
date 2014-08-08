/*
 * tga_writer.h
 *
 *  Created on: Apr 27, 2013
 *      Author: Sachith Dunatunga
 */

#ifndef __TGA_WRITER_H__
#define __TGA_WRITER_H__

typedef struct s_tgaheader {
	char  idlength;
	char  colourmaptype;
	char  datatypecode;
	short int colourmaporigin;
	short int colourmaplength;
	char  colourmapdepth;
	short int x_origin;
	short int y_origin;
	short width;
	short height;
	char  bitsperpixel;
	char  imagedescriptor;
	char  idfield[255];
} tgaheader_t;

void create_tga_header(tgaheader_t *th, int width, int height, int bitdepth);
void write_tga_header(FILE *fp, tgaheader_t *th);
void write_tga_image(FILE *fp, tgaheader_t *th, char *tga_pixeldata);

#endif /* __TGA_WRITER_H__ */
