/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_inoutput.h : input and output of matrices.               */
/*     									   */
/* Date : July 8, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Author nor the Institution (National Institute of Standards
// and Technology) make any representations about the suitability of this 
// software for any purpose. This software is provided "as is" without 
// expressed or implied warranty.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef __GMM_INOUTPUT_H
#define __GMM_INOUTPUT_H

namespace gmm {

  /*************************************************************************/
  /*                                                                       */
  /*  Functions to read and write Harwell Boeing format.                   */
  /*                                                                       */
  /*************************************************************************/

  inline void IOHBTerminate(const char *a) { DAL_THROW(failure_error, a); }

  inline char* substr(const char* S, int pos, int len) {
    int i;
    char *SubS;
    if ( pos+len <= int(strlen(S))) {
      SubS = (char *)malloc(len+1);
      if ( SubS == 0 ) IOHBTerminate("Insufficient memory for SubS.");
      for (i=0;i<len;i++) SubS[i] = S[pos+i];
      SubS[len] = (char) NULL;
    } else {
      SubS = NULL;
    }
    return SubS;
  }

  inline void upcase(char *t)
  { if (t) while (*t) { *t = toupper(*t); ++t; } }

  inline int ParseIfmt(char* fmt, int* perline, int* width) {
    /*************************************************/
    /*  Parse an *integer* format field to determine */
    /*  width and number of elements per line.       */
    /*************************************************/
    char *tmp;
    if (fmt == NULL ) {
      *perline = 0; *width = 0; return 0;
    }
    tmp = strchr(fmt,'(');
    tmp = substr(fmt,tmp - fmt + 1, strchr(fmt,'I') - tmp - 1);
    *perline = atoi(tmp);
    tmp = strchr(fmt,'I');
    tmp = substr(fmt,tmp - fmt + 1, strchr(fmt,')') - tmp - 1);
    return *width = atoi(tmp);
  }
  
  inline int ParseRfmt(char* fmt, int* perline, int* width,
		       int* prec, int* flag) {
    /*************************************************/
    /*  Parse a *real* format field to determine     */
    /*  width and number of elements per line.       */
    /*  Also sets flag indicating 'E' 'F' 'P' or 'D' */
    /*  format.                                      */
    /*************************************************/
    char* tmp;
    char* tmp2;
    char* tmp3;
    int len;
    
    if (fmt == NULL ) {
      *perline = 0; 
      *width = 0; 
      flag = NULL;  
      return 0;
    }
    
    if (strchr(fmt,'(') != NULL)  fmt = strchr(fmt,'(');
    if (strchr(fmt,')') != NULL)  {
      tmp2 = strchr(fmt,')');
      while ( strchr(tmp2+1,')') != NULL ) {
	tmp2 = strchr(tmp2+1,')');
      }
      // *(tmp2+1) = 0;
    }
    if (strchr(fmt,'P') != NULL)  /* Remove any scaling factor, which */
      {                             /* affects output only, not input */
	if (strchr(fmt,'(') != NULL) {
	  tmp = strchr(fmt,'P');
	  if ( *(++tmp) == ',' ) tmp++;
	  tmp3 = strchr(fmt,'(')+1;
	  len = tmp-tmp3;
	  tmp2 = tmp3;
	  while ( *(tmp2+len) != 0 ) {
	    *tmp2=*(tmp2+len);
	    tmp2++; 
	  }
	  // *(strchr(fmt,')')+1) = 0;
	}
      }
    if (strchr(fmt,'E') != NULL) { 
      *flag = 'E';
    } else if (strchr(fmt,'D') != NULL) { 
      *flag = 'D';
    } else if (strchr(fmt,'F') != NULL) { 
      *flag = 'F';
    } else {
      fprintf(stderr,"Real format %s in H/B file not supported.\n",fmt);
      return 0;
    }
    tmp = strchr(fmt,'(');
    tmp = substr(fmt,tmp - fmt + 1, strchr(fmt,*flag) - tmp - 1);
    *perline = atoi(tmp);
    tmp = strchr(fmt,*flag);
    if ( strchr(fmt,'.') ) {
      *prec = atoi(substr( fmt, strchr(fmt,'.') - fmt + 1, strchr(fmt,')') 
			   - strchr(fmt,'.')-1) );
      tmp = substr(fmt,tmp - fmt + 1, strchr(fmt,'.') - tmp - 1);
    } else {
      tmp = substr(fmt,tmp - fmt + 1, strchr(fmt,')') - tmp - 1);
    }
    return *width = atoi(tmp);
  }
  
  
  inline int readHB_header(FILE* in_file, char* Title,
			   char* Key, char* Type, 
			   int* Nrow, int* Ncol, int* Nnzero,
			   int* Nrhs, char* Ptrfmt,
			   char* Indfmt, char* Valfmt,
			   char* Rhsfmt, int* Ptrcrd,
			   int* Indcrd, int* Valcrd,
			   int* Rhscrd, char *Rhstype) {
    /*************************************************************************/
    /*  Read header information from the named H/B file...                   */
    /*************************************************************************/
    int Totcrd,Neltvl,Nrhsix;
    char line[BUFSIZ];
    
    /*  First line:   */
    fgets(line, BUFSIZ, in_file);
    if ( sscanf(line,"%*s") < 0 ) 
      IOHBTerminate("iohb.c: Null (or blank) first line of HB file.\n");
    (void) sscanf(line, "%72c%8[^\n]", Title, Key);
    *(Key+8) = (char) NULL;
    *(Title+72) = (char) NULL;
    
    /*  Second line:  */
    fgets(line, BUFSIZ, in_file);
    if ( sscanf(line,"%*s") < 0 ) 
      IOHBTerminate("iohb.c: Null (or blank) second line of HB file.\n");
    if ( sscanf(line,"%i",&Totcrd) != 1) Totcrd = 0;
    if ( sscanf(line,"%*i%i",Ptrcrd) != 1) *Ptrcrd = 0;
    if ( sscanf(line,"%*i%*i%i",Indcrd) != 1) *Indcrd = 0;
    if ( sscanf(line,"%*i%*i%*i%i",Valcrd) != 1) *Valcrd = 0;
    if ( sscanf(line,"%*i%*i%*i%*i%i",Rhscrd) != 1) *Rhscrd = 0;
    
    /*  Third line:   */
    fgets(line, BUFSIZ, in_file);
    if ( sscanf(line,"%*s") < 0 ) 
      IOHBTerminate("Null (or blank) third line of HB file.\n");
    if ( sscanf(line, "%3c", Type) != 1) 
      IOHBTerminate("Invalid Type info, line 3 of Harwell-Boeing file.\n");
    upcase(Type);
    if ( sscanf(line,"%*3c%i",Nrow) != 1) *Nrow = 0 ;
    if ( sscanf(line,"%*3c%*i%i",Ncol) != 1) *Ncol = 0 ;
    if ( sscanf(line,"%*3c%*i%*i%i",Nnzero) != 1) *Nnzero = 0 ;
    if ( sscanf(line,"%*3c%*i%*i%*i%i",&Neltvl) != 1) Neltvl = 0 ;
    
    /*  Fourth line:  */
    fgets(line, BUFSIZ, in_file);
    if ( sscanf(line,"%*s") < 0 ) 
      IOHBTerminate("Null (or blank) fourth line of HB file.\n");
    if ( sscanf(line, "%16c",Ptrfmt) != 1)
      IOHBTerminate("Invalid format info, line 4 of Harwell-Boeing file.\n"); 
    if ( sscanf(line, "%*16c%16c",Indfmt) != 1)
      IOHBTerminate("Invalid format info, line 4 of Harwell-Boeing file.\n"); 
    if ( sscanf(line, "%*16c%*16c%20c",Valfmt) != 1) 
      IOHBTerminate("Invalid format info, line 4 of Harwell-Boeing file.\n"); 
    sscanf(line, "%*16c%*16c%*20c%20c",Rhsfmt);
    *(Ptrfmt+16) = (char) NULL;
    *(Indfmt+16) = (char) NULL;
    *(Valfmt+20) = (char) NULL;
    *(Rhsfmt+20) = (char) NULL;
    
    /*  (Optional) Fifth line: */
    if (*Rhscrd != 0 )
      { 
	fgets(line, BUFSIZ, in_file);
	if ( sscanf(line,"%*s") < 0 ) 
	  IOHBTerminate("Null (or blank) fifth line of HB file.\n");
	if ( sscanf(line, "%3c", Rhstype) != 1) 
	  IOHBTerminate("Invalid RHS type information, line 5 of"
			" Harwell-Boeing file.\n");
	if ( sscanf(line, "%*3c%i", Nrhs) != 1) *Nrhs = 0;
	if ( sscanf(line, "%*3c%*i%i", &Nrhsix) != 1) Nrhsix = 0;
      }
    return 1;
  }
  
  inline int readHB_info(const char* filename, int* M, 
		  int* N, int* nz, char** Type, int* Nrhs) {
    /***********************************************************************/
    /*  The readHB_info function opens and reads the header information    */
    /*  form the specified Harwell-Boeing file, and reports back the       */
    /*  number of rows and columns in the stored matrix (M and N), the     */
    /*  number of nonzeros in the matrix (nz), and the number of           */
    /*  right-hand-sides stored along with the matrix (Nrhs).              */
    /*                                                                     */
    /*  For a description of the Harwell Boeing standard, see:             */
    /*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989         */
    /*                                                                     */
    /*    ----------                                                       */
    /*    **CAVEAT**                                                       */
    /*    ----------                                                       */
    /*  **  If the input file does not adhere to the H/B format, the  **   */
    /*  **             results will be unpredictable.                 **   */
    /*                                                                     */
    /***********************************************************************/
    FILE *in_file;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd; 
    int Nrow, Ncol, Nnzero;
    char* mat_type;
    char Title[73], Key[9], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];

    mat_type = *Type;
    if ( mat_type == NULL ) IOHBTerminate("Insufficient memory for mat_typen");
    
    if ( (in_file = fopen( filename, "r")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
    }
    
    readHB_header(in_file, Title, Key, mat_type, &Nrow, &Ncol, &Nnzero, Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt, 
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
    fclose(in_file);
    *Type = mat_type;
    *(*Type+3) = (char) NULL;
    *M    = Nrow;
    *N    = Ncol;
    *nz   = Nnzero;
    if (Rhscrd == 0) {*Nrhs = 0;}
    return 1;
    
  }
  

  
  inline int readHB_mat_double(const char* filename,
			       size_type colptr[], size_type rowind[], 
			       double val[]) {
    /************************************************************************/
    /*  This function opens and reads the specified file, interpreting its  */
    /*  contents as a sparse matrix stored in the Harwell/Boeing standard   */
    /*  format and creating compressed column storage scheme vectors to hold*/
    /*  the index and nonzero value information.                            */
    /*                                                                      */
    /*    ----------                                                        */
    /*    **CAVEAT**                                                        */
    /*    ----------                                                        */
    /*  Parsing real formats from Fortran is tricky, and this file reader   */
    /*  does not claim to be foolproof.   It has been tested for cases when */
    /*  the real values are printed consistently and evenly spaced on each  */
    /*  line, with Fixed (F), and Exponential (E or D) formats.             */
    /*                                                                      */
    /*  **  If the input file does not adhere to the H/B format, the  **    */
    /*  **             results will be unpredictable.                 **    */
    /*                                                                      */
    /************************************************************************/
    FILE *in_file;
    int i,j,ind,col,offset,count,last,Nrhs;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero, Nentries;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Valperline, Valwidth, Valprec;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    char* ThisElement;
    char Title[73], Key[8], Type[4], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    char line[BUFSIZ];
    
    if ( (in_file = fopen( filename, "r")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
    }
    
    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
    
    /*  Parse the array input formats from Line 3 of HB file  */
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
      ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
    }
    
    /*  Read column pointer array:   */
    
    offset = 0;         /* if base 0 storage is declared (via macro def),  */
                        /* then storage entries are offset by 1            */
    
    ThisElement = (char *) malloc(Ptrwidth+1);
    if ( ThisElement == NULL )
      IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Ptrwidth) = (char) NULL;
    count=0;
    for (i=0;i<Ptrcrd;i++) {
      fgets(line, BUFSIZ, in_file);
      if ( sscanf(line,"%*s") < 0 ) 
	IOHBTerminate("Null (or blank) line in pointer data region.\n");
      col =  0;
      for (ind = 0;ind<Ptrperline;ind++) {
	if (count > Ncol) break;
	strncpy(ThisElement,line+col,Ptrwidth);
	/* ThisElement = substr(line,col,Ptrwidth); */
	colptr[count] = atoi(ThisElement)-offset;
	count++; col += Ptrwidth;
      }
    }
    free(ThisElement);
    
    /*  Read row index array:  */
    
    ThisElement = (char *) malloc(Indwidth+1);
    if ( ThisElement == NULL )
      IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Indwidth) = (char) NULL;
    count = 0;
    for (i=0;i<Indcrd;i++) {
      fgets(line, BUFSIZ, in_file);
      if ( sscanf(line,"%*s") < 0 ) 
	IOHBTerminate("Null (or blank) line in index data region.\n");
      col =  0;
      for (ind = 0;ind<Indperline;ind++) {
	if (count == Nnzero) break;
	strncpy(ThisElement,line+col,Indwidth);
	/*        ThisElement = substr(line,col,Indwidth); */
	rowind[count] = atoi(ThisElement)-offset;
	count++; col += Indwidth;
      }
    }
    free(ThisElement);
    
    /*  Read array of values:  */
    
    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
      
      if ( Type[0] == 'C' ) Nentries = 2*Nnzero;
      else Nentries = Nnzero;
      
      ThisElement = (char *) malloc(Valwidth+1);
      if ( ThisElement == NULL )
	IOHBTerminate("Insufficient memory for ThisElement.");
      *(ThisElement+Valwidth) = (char) NULL;
      count = 0;
      for (i=0;i<Valcrd;i++) {
	fgets(line, BUFSIZ, in_file);
	if ( sscanf(line,"%*s") < 0 ) 
	  IOHBTerminate("Null (or blank) line in value data region.\n");
	if (Valflag == 'D')  {
          while( strchr(line,'D') ) *strchr(line,'D') = 'E';
	  /*           *strchr(Valfmt,'D') = 'E'; */
	}
	col =  0;
	for (ind = 0;ind<Valperline;ind++) {
	  if (count == Nentries) break;
	  strncpy(ThisElement,line+col,Valwidth);
          /*ThisElement = substr(line,col,Valwidth);*/
          if ( Valflag != 'F' && strchr(ThisElement,'E') == NULL ) { 
	    /* insert a char prefix for exp */
	    last = strlen(ThisElement);
	    for (j=last+1;j>=0;j--) {
	      ThisElement[j] = ThisElement[j-1];
	      if ( ThisElement[j] == '+' || ThisElement[j] == '-' ) {
		ThisElement[j-1] = Valflag;                    
		break;
	      }
	    }
          }
          val[count] = atof(ThisElement);
          count++; col += Valwidth;
	}
      }
      free(ThisElement);
    }
    
    fclose(in_file);
    return 1;
  }


  inline int writeHB_mat_double(const char* filename,
				int M, int N, int nz,
				const size_type colptr[],
				const size_type rowind[], 
				const double val[], int Nrhs,
				const double rhs[], 
				const double guess[],
				const double exact[],
				const char* Title,
				const char* Key,
				const char* Type, 
			        char* Ptrfmt, char* Indfmt,
			        char* Valfmt, char* Rhsfmt,
				const char* Rhstype, int shift) {
    /************************************************************************/
    /*  The writeHB function opens the named file and writes the specified  */
    /*  matrix and optional right-hand-side(s) to that file in              */
    /*  Harwell-Boeing format.                                              */
    /*                                                                      */
    /*  For a description of the Harwell Boeing standard, see:              */
    /*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989          */
    /*                                                                      */
    /************************************************************************/
    FILE *out_file;
    int i,j,entry,offset,acount,linemod;
    int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
    int nvalentries, nrhsentries;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Rhsperline, Rhswidth, Rhsprec;
    int Rhsflag;
    int Valperline, Valwidth, Valprec;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    char pformat[16],iformat[16],vformat[19],rformat[19];
    
    if ( Type[0] == 'C' ) {
      nvalentries = 2*nz;
      nrhsentries = 2*M;
    } else {
      nvalentries = nz;
      nrhsentries = M;
    }
    
    if ( filename != NULL ) {
      if ( (out_file = fopen( filename, "w")) == NULL ) {
	fprintf(stderr,"Error: Cannot open file: %s\n",filename);
	return 0;
      }
    } else out_file = stdout;
    
    if ( Ptrfmt == NULL ) Ptrfmt = const_cast<char *>("(8I10)");
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    sprintf(pformat,"%%%dd",Ptrwidth);
    ptrcrd = (N+1)/Ptrperline;
    if ( (N+1)%Ptrperline != 0) ptrcrd++;
    
    if ( Indfmt == NULL ) Indfmt =  Ptrfmt;
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    sprintf(iformat,"%%%dd",Indwidth);
    indcrd = nz/Indperline;
    if ( nz%Indperline != 0) indcrd++;
    
    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
      if ( Valfmt == NULL ) Valfmt = const_cast<char *>("(4E20.13)");
      ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
      if (Valflag == 'D') *strchr(Valfmt,'D') = 'E';
      if (Valflag == 'F')
	sprintf(vformat,"%% %d.%df",Valwidth,Valprec);
      else
	sprintf(vformat,"%% %d.%dE",Valwidth,Valprec);
      valcrd = nvalentries/Valperline;
      if ( nvalentries%Valperline != 0) valcrd++;
    } else valcrd = 0;
    
    if ( Nrhs > 0 ) {
      if ( Rhsfmt == NULL ) Rhsfmt = Valfmt;
      ParseRfmt(Rhsfmt,&Rhsperline,&Rhswidth,&Rhsprec, &Rhsflag);
      if (Rhsflag == 'F')
	sprintf(rformat,"%% %d.%df",Rhswidth,Rhsprec);
      else
	sprintf(rformat,"%% %d.%dE",Rhswidth,Rhsprec);
      if (Rhsflag == 'D') *strchr(Rhsfmt,'D') = 'E';
      rhscrd = nrhsentries/Rhsperline; 
      if ( nrhsentries%Rhsperline != 0) rhscrd++;
      if ( Rhstype[1] == 'G' ) rhscrd+=rhscrd;
      if ( Rhstype[2] == 'X' ) rhscrd+=rhscrd;
      rhscrd*=Nrhs;
    } else rhscrd = 0;
    
    totcrd = 4+ptrcrd+indcrd+valcrd+rhscrd;
    
    
    /*  Print header information:  */
    
    fprintf(out_file,"%-72s%-8s\n%14d%14d%14d%14d%14d\n",Title, Key, totcrd,
            ptrcrd, indcrd, valcrd, rhscrd);
    fprintf(out_file,"%3s%11s%14d%14d%14d\n",Type,"          ", M, N, nz);
    fprintf(out_file,"%-16s%-16s%-20s", Ptrfmt, Indfmt, Valfmt);
    if ( Nrhs != 0 ) {
      /*    Print Rhsfmt on fourth line and                                 */
      /*      optional fifth header line for auxillary vector information:  */
      fprintf(out_file,"%-20s\n%-14s%d\n",Rhsfmt,Rhstype,Nrhs);
    } else fprintf(out_file,"\n");
    
    offset = 1 - shift;  /* if base 0 storage is declared (via macro def), */
                         /* then storage entries are offset by 1           */
    
    /*  Print column pointers:   */
    for (i=0;i<N+1;i++) {
      entry = colptr[i]+offset;
      fprintf(out_file,pformat,entry);
      if ( (i+1)%Ptrperline == 0 ) fprintf(out_file,"\n");
    }
    
    if ( (N+1) % Ptrperline != 0 ) fprintf(out_file,"\n");
    
    /*  Print row indices:       */
    for (i=0;i<nz;i++) {
      entry = rowind[i]+offset;
      fprintf(out_file,iformat,entry);
      if ( (i+1)%Indperline == 0 ) fprintf(out_file,"\n");
    }
    
    if ( nz % Indperline != 0 ) fprintf(out_file,"\n");
    
    /*  Print values:            */
    
    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
      
      for (i=0;i<nvalentries;i++) {
	fprintf(out_file,vformat,val[i]);
	if ( (i+1)%Valperline == 0 ) fprintf(out_file,"\n");
      }
      
      if ( nvalentries % Valperline != 0 ) fprintf(out_file,"\n");
      
      /*  If available,  print right hand sides, 
	  guess vectors and exact solution vectors:  */
      acount = 1;
      linemod = 0;
      if ( Nrhs > 0 ) {
	for (i=0;i<Nrhs;i++) {
          for ( j=0;j<nrhsentries;j++ ) {
            fprintf(out_file,rformat,rhs[j]);
            if ( acount++%Rhsperline == linemod ) fprintf(out_file,"\n");
          }
          if ( acount%Rhsperline != linemod ) {
            fprintf(out_file,"\n");
            linemod = (acount-1)%Rhsperline;
          }
          rhs += nrhsentries;
          if ( Rhstype[1] == 'G' ) {
            for ( j=0;j<nrhsentries;j++ ) {
              fprintf(out_file,rformat,guess[j]);
              if ( acount++%Rhsperline == linemod ) fprintf(out_file,"\n");
            }
            if ( acount%Rhsperline != linemod ) {
              fprintf(out_file,"\n");
              linemod = (acount-1)%Rhsperline;
            }
            guess += nrhsentries;
          }
          if ( Rhstype[2] == 'X' ) {
            for ( j=0;j<nrhsentries;j++ ) {
              fprintf(out_file,rformat,exact[j]);
              if ( acount++%Rhsperline == linemod ) fprintf(out_file,"\n");
            }
            if ( acount%Rhsperline != linemod ) {
              fprintf(out_file,"\n");
              linemod = (acount-1)%Rhsperline;
            }
            exact += nrhsentries;
          }
	}
      }
      
    }
    
    if ( fclose(out_file) != 0){
      fprintf(stderr,"Error closing file in writeHB_mat_double().\n");
      return 0;
    } else return 1;
    
  }

  // not securized, to be used with "double" or "std::complex<double>"

  template <class T, int shift> void
  Harwell_Boeing_save(const char *filename, const csc_matrix<T, shift>& A) {
    static const char *type1 = "RRA";  // real not squared
    static const char *type2 = "RUA";  // real squared
    static const char *type3 = "CRA";  // complex not squared
    static const char *type4 = "CUA";  // complex squared
    const char *t = 0;

    if (sizeof(T) == sizeof(std::complex<double>))
      if (mat_nrows(A) == mat_ncols(A)) t = type4; else t = type3;
    else
      if (mat_nrows(A) == mat_ncols(A)) t = type2; else t = type1;
    
    writeHB_mat_double(filename, mat_nrows(A), mat_ncols(A),
		       A.jc[mat_ncols(A)], A.jc, A.ir,
		       (double *)A.pr,
		       0, 0, 0, 0, "GETFEM++ CSC MATRIX", "CSCMAT",
		       t, 0, 0, 0, 0, "F", shift);
  }

  template <class T, int shift> void
  Harwell_Boeing_load(const char *filename, csc_matrix<T, shift>& A) {
    int M, N, nonzeros;
    int Nrhs;
    char* Type;
    Type = new char[4];
    readHB_info(filename, &M, &N, &nonzeros, &Type, &Nrhs);

    if (A.pr) { delete[] A.pr; delete[] A.ir; delete[] A.jc; }
    A.nc = N; A.nr = M;
    A.jc = new size_type[N+1];
    A.ir = new size_type[nonzeros];
    if (!(A.jc) || !(A.ir))
      DAL_THROW(failure_error, "Insufficient memory for rowind.\n");
    
    if (Type[0] == 'C') {
      A.pr = new double[2*nonzeros];
      if (!(A.pr)) DAL_THROW(failure_error, "Insufficient memory for val.\n");
      
    } else {
      if (Type[0] != 'P') { 
	A.pr = new double[nonzeros];
	if (!(A.pr))
	  DAL_THROW(failure_error, "Insufficient memory for val.\n");
      }
    }
    readHB_mat_double(filename, A.jc, A.ir, (double *)(A.pr));
    for (size_type i = 0; i <= N; ++i)
      { A.jc[i] += shift; A.jc[i] -= 1; }
    for (size_type i = 0; i < nonzeros; ++i)
      { A.ir[i] += shift; A.ir[i] -= 1; }
    delete [] Type;
  }



}


#endif //  __GMM_INOUTPUT_H
