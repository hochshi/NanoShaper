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

#include "../globals.h"
#include "ply.h"

// Swap endian-ness for four-byte elements

inline void endian_swap_long(unsigned char *p)
{
 unsigned char b0,b1,b2,b3;

 b0 = *p; b1 = *(p+1); b2 = *(p+2); b3 = *(p+3);
 *p = b3; *(p+1) = b2; *(p+2) = b1; *(p+3) = b0;
}


// Read one line (max 1024 chars) and exit if EOF

char *readLineFromFile(FILE *in, bool exit_on_eof)
{
#define MAX_READLINE_CHARS	1024
 static char line[MAX_READLINE_CHARS];
 int i=0;
 char c;

 while ((c = fgetc(in)) != '\n' && i<(MAX_READLINE_CHARS-1))
   if (c==EOF)
   {
    if (exit_on_eof) cout << endl << ERR <<("\nUnexpected end of file!\n");
    else return NULL;
   }
   else if (c != '\r') line[i++] = c;
 line[i] = '\0';

 if (i==MAX_READLINE_CHARS-1)
  cout << endl << WARN <<("readLineFromFile: Line is too long. Truncated !\n");

 return line;
}


// Looks for a keyword 'kw' in an ASCII file referenced through 'fp'.
// The file pointer is set to the byte right after the first keyword matched.
// Return 1 on success (keyword match), 0 otherwise.


bool seek_keyword(FILE *fp, const char *kw)
{
 static char s[256];
 s[0]='\0';
 do 
 {
   int check = fscanf(fp,"%255s",s); 
   if (check<=0)
	   return 0;
 } while (strcmp(s,kw) && !feof(fp));
 if (feof(fp)) return 0;
 return 1;
}


inline void skipCommentAndBlankLines(FILE *fp)
{
 long pos0;
 char *line, s[2];
 do {pos0 = ftell(fp); line = readLineFromFile(fp);} while (line[0] == '#' || line[0] == '\0' || !sscanf(line,"%1s",s));
 fseek(fp, pos0, SEEK_SET);
}

////////////////////// PLY LOADER //////////////////////////////////////////////

int ply_parseElements(FILE *in, const char *elname)
{
 char c, keyword[64];
 int num;
 // skip comments
 if (!fscanf(in,"%64s ",keyword)) cout << endl << ERR <<("Unexpected token or end of file!\n");
 while (!strcmp(keyword,"comment") || !strcmp(keyword,"obj_info"))
 {
  while ((c = fgetc(in)) != '\n') if (c==EOF) cout << endl << ERR <<("\nUnexpected end of file!\n");
  if (!fscanf(in,"%64s ",keyword)) cout << endl << ERR <<("Unexpected token or end of file!\n");
 }
 if (strcmp(keyword,"element")) cout << endl << ERR <<("element definition expected!\n");
 if (!fscanf(in,"%64s ",keyword)) cout << endl << ERR <<("Unexpected token or end of file!\n");
 if (strcmp(keyword,elname)) cout << endl << ERR <<("Sorry. Element type '%s' is not supported!\n",keyword);
 if (!fscanf(in,"%d\n",&num)) cout << endl << ERR <<("Unexpected token or end of file!\n");
 if (num <= 0) cout << endl << ERR <<("Unexpected empty element list!\n");

 return num;
}

void ply_checkVertexProperties(FILE *in)
{
 char keyword[64], dtype[64], dval[64];
 if (fscanf(in,"%64s %64s %64s\n",keyword,dtype,dval) < 3) cout << endl << ERR <<("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) cout << endl << ERR <<("property definition expected!\n");
 if (strcmp(dtype,"float") && strcmp(dtype,"float32")) cout << endl << ERR <<("float property expected!\n");
 if (strcmp(dval,"x")) cout << endl << ERR <<("'x' float property expected!\n");
 if (fscanf(in,"%64s %64s %64s\n",keyword,dtype,dval) < 3) cout << endl << ERR <<("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) cout << endl << ERR <<("property definition expected!\n");
 if (strcmp(dtype,"float") && strcmp(dtype,"float32")) cout << endl << ERR <<("float property expected!\n");
 if (strcmp(dval,"y")) cout << endl << ERR <<("'y' float property expected!\n");
 if (fscanf(in,"%64s %64s %64s\n",keyword,dtype,dval) < 3) cout << endl << ERR <<("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) cout << endl << ERR <<("property definition expected!\n");
 if (strcmp(dtype,"float") && strcmp(dtype,"float32")) cout << endl << ERR <<("float property expected!\n");
 if (strcmp(dval,"z")) cout << endl << ERR <<("'z' float property expected!\n");
}

int ply_getOverhead(FILE *in, int format, const char *element)
{
 char keyword[64], ptype[64], pname[64];
 int oh = 0;
 long pos = ftell(in);
 char *rline = readLineFromFile(in);
 if (!sscanf(rline,"%64s ",keyword)) cout << endl << ERR <<("Unexpected token or end of file!\n");
 while (!strcmp(keyword, "property"))
 {
  if (sscanf(rline,"%64s %64s %64s",keyword,ptype,pname) < 3) cout << endl << ERR <<("Unexpected token or end of file!\n");
  if (!strcmp(element,"vertex") && !strcmp(pname,"x")) break;
  else if (!strcmp(element,"face") && !strcmp(ptype,"list")) break;
  pos = ftell(in);
  if (!strcmp(ptype, "char") || !strcmp(ptype, "uchar")) oh += (format)?(1):1;
  else if (!strcmp(ptype, "short") || !strcmp(ptype, "ushort")) oh += (format)?(2):1;
  else if (!strcmp(ptype, "int") || !strcmp(ptype, "uint") || 
           !strcmp(ptype, "float") || !strcmp(ptype,"float32")) oh += (format)?(4):1;
  else if (!strcmp(ptype, "double")) oh += (format)?(8):1;
  else if (!strcmp(ptype, "list")) cout << endl << ERR <<("list properties other than face indices are not supported!\n");
  else cout << endl << ERR <<("Unrecognized property type!\n");
  if (!sscanf(readLineFromFile(in),"%64s ",keyword)) cout << endl << ERR <<("Unexpected token or end of file!\n");
 }
 fseek(in, pos, SEEK_SET);

 return oh;
}

void ply_checkFaceProperties(FILE *in)
{
 char keyword[64], ltype[64], uctype[64], dtype[64], dval[64];
 if (fscanf(in,"%64s %64s %64s %64s %64s\n",keyword,ltype,uctype,dtype,dval) < 5) cout << endl << ERR <<("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) cout << endl << ERR <<("property definition expected!\n");
 if (strcmp(ltype,"list")) cout << endl << ERR <<("list property expected!\n");
 if (strcmp(uctype,"uchar") && strcmp(uctype,"uint8")) cout << endl << ERR <<("uchar property expected!\n");
 if (strcmp(dtype,"int") && strcmp(dtype,"int32")) cout << endl << ERR <<("int property expected!\n");
 if (strcmp(dval,"vertex_indices")) cout << endl << ERR <<("vertex_indices property expected!\n");
}

void ply_readOverhead(FILE *in, int format, int oh)
{
 int i;
 static char token[1024];
 if (format == PLY_FORMAT_ASCII) 
	 for (i=0; i<oh; i++) 
	 {
		 int check = fscanf(in, "%s", token);
		 if (check<=0)
		 {
			cout << endl << ERR << "Error in reading ply file";
			cout << endl;
			exit(-1);
		 }
	 }
 else for (i=0; i<oh; i++) fgetc(in);
}


int ply_readVCoords(FILE *in, int format, int ph, int oh, float *x, float *y, float *z)
{
 float vc[3];

 ply_readOverhead(in, format, ph);

 if (format == PLY_FORMAT_ASCII)
 {
  if (fscanf(in,"%f %f %f", x, y, z) < 3) cout << endl << ERR <<("Unexpected token or end of file!\n"); 
 }
 else
 {
  if (fread(vc, 4, 3, in) < 3) cout << endl << ERR <<("Unexpected end of file!\n");
  *x = vc[0]; *y = vc[1]; *z = vc[2];

  if (format == PLY_FORMAT_BIN_B)
  {
   endian_swap_long((unsigned char *)(x));
   endian_swap_long((unsigned char *)(y));
   endian_swap_long((unsigned char *)(z));
  }
 }

 ply_readOverhead(in, format, oh);

 return 1;
}

int ply_readFIndices(FILE *in, int format, int ph, int *nv, int *x, int *y, int *z)
{
 unsigned char nvs;
 int vc[3];

 ply_readOverhead(in, format, ph);

 if (format == PLY_FORMAT_ASCII) 
 {
	 int check = fscanf(in,"%d %d %d %d", nv, x, y, z); 
	 if (check<=0)
	 {
		cout << endl << ERR << "Error in reading ply file!";
		cout << endl;
		exit(-1);
	 }
	 return 1;
 }

 size_t check = fread(&nvs, 1, 1, in);
 if (check<=(size_t)0)
 {
 	cout << endl << ERR << "Error in reading ply file!";
	cout << endl;
	exit(-1);
 }
 *nv = (int)nvs;
 check = fread(vc, 4, 3, in);
 if (check<=(size_t)0)
 {
 	cout << endl << ERR << "Error in reading ply file!";
	cout << endl;
	exit(-1);
 }
 *x = vc[0]; *y = vc[1]; *z = vc[2];

 if (format == PLY_FORMAT_BIN_B)
 {
  endian_swap_long((unsigned char *)(x));
  endian_swap_long((unsigned char *)(y));
  endian_swap_long((unsigned char *)(z));
 }

 return 1;
}

int ply_readAnotherFIndex(FILE *in, int format, int *x)
{
 if (format == PLY_FORMAT_ASCII) return (fscanf(in,"%d", x));

 if (fread(x, 4, 1, in) < 1) cout << endl << ERR <<("Unexpected end of file!\n");

 if (format == PLY_FORMAT_BIN_B) endian_swap_long((unsigned char *)(x));

 return 1;
}

