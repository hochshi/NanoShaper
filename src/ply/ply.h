/****************************************************************************
* JMeshLib                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#define PLY_FORMAT_ASCII	0
#define PLY_FORMAT_BIN_L	1
#define PLY_FORMAT_BIN_B	2

#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW

inline void endian_swap_long(unsigned char *p);
char *readLineFromFile(FILE *in, bool exit_on_eof = 1);
bool seek_keyword(FILE *fp, const char *kw);
inline void skipCommentAndBlankLines(FILE *fp);
int ply_parseElements(FILE *in, const char *elname);
void ply_checkVertexProperties(FILE *in);
int ply_getOverhead(FILE *in, int format, const char *element);
void ply_checkFaceProperties(FILE *in);
void ply_readOverhead(FILE *in, int format, int oh);
int ply_readVCoords(FILE *in, int format, int ph, int oh, float *x, float *y, float *z);
int ply_readFIndices(FILE *in, int format, int ph, int *nv, int *x, int *y, int *z);
int ply_readAnotherFIndex(FILE *in, int format, int *x);