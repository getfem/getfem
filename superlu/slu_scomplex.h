
/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
*/
/*
  Copyright (c) 1997 by Xerox Corporation.  All rights reserved.
 
  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 
  Permission is hereby granted to use or copy this program for any
  purpose, provided the above notices are retained on all copies.
  Permission to modify the code and to distribute modified code is
  granted, provided the above notices are retained, and a notice that
  the code was modified is included with the above copyright notice.
*/
#ifndef __SUPERLU_SCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_SCOMPLEX

/* 
 * This header file is to be included in source files c*.c
 */
#ifndef SCOMPLEX_INCLUDE
#define SCOMPLEX_INCLUDE

typedef struct { float r, i; } complex;


/* Macro definitions */

/* Complex Addition c = a + b */
#define c_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/* Complex Subtraction c = a - b */
#define c_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/* Complex-Double Multiplication */
#define cs_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/* Complex-Complex Multiplication */
#define cc_mult(c, a, b) { \
	float cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

#define cc_conj(a, b) { \
        (a)->r = (b)->r; \
        (a)->i = -((b)->i); \
    }

/* Complex equality testing */
#define c_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in scomplex.c */
void c_div(complex *, complex *, complex *);
double c_abs(complex *);     /* exact */
double c_abs1(complex *);    /* approximate */
void c_exp(complex *, complex *);
void r_cnjg(complex *, complex *);
double r_imag(complex *);


#ifdef __cplusplus
  }
#endif

#endif

#endif  /* __SUPERLU_SCOMPLEX */
