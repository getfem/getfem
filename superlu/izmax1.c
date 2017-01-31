#include <math.h>
#include "slu_Cnames.h"
#include "slu_dcomplex.h"

int
izmax1_(int *n, doublecomplex *cx, int *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Copyright (c) 1992-2013 The University of Tennessee and The University
                        of Tennessee Research Foundation.  All rights
                        reserved.
       Copyright (c) 2000-2013 The University of California Berkeley. All
                        rights reserved.
       Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
                        reserved.

       Redistribution and use in source and binary forms, with or without
       modification, are permitted provided that the following conditions are
       met:

       - Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.

       - Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer listed
         in this license in the documentation and/or other materials
         provided with the distribution.

       - Neither the name of the copyright holders nor the names of its
         contributors may be used to endorse or promote products derived from
         this software without specific prior written permission.

       The copyright holders provide no reassurances that the source code
       provided does not infringe any patent, copyright, or any other
       intellectual property rights of third parties.  The copyright holders
       disclaim any liability to any recipient for claims brought against
       recipient by any third party for infringement of that parties
       intellectual property rights.

       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
       "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
       LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
       A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
       OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
       SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
       LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
       THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
       (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
       OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


    Purpose   
    =======   

    IZMAX1 finds the index of the element whose real part has maximum   
    absolute value.   

    Based on IZAMAX from Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with ZLACON.   

    Arguments   
    =========   

    N       (input) INT   
            The number of elements in the vector CX.   

    CX      (input) COMPLEX*16 array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INT   
            The spacing between successive values of CX.  INCX >= 1.   

   ===================================================================== 
*/  

    /* System generated locals */
    int ret_val, i__1, i__2;
    double d__1;
    
    /* Local variables */
    double smax;
    int i, ix;

#define CX(I) cx[(I)-1]

    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L30;
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    smax = (d__1 = CX(1).r, fabs(d__1));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = ix;
	if ((d__1 = CX(ix).r, fabs(d__1)) <= smax) {
	    goto L10;
	}
	ret_val = i;
	i__2 = ix;
	smax = (d__1 = CX(ix).r, fabs(d__1));
L10:
	ix += *incx;
/* L20: */
    }
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

L30:
    smax = (d__1 = CX(1).r, fabs(d__1));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = i;
	if ((d__1 = CX(i).r, fabs(d__1)) <= smax) {
	    goto L40;
	}
	ret_val = i;
	i__2 = i;
	smax = (d__1 = CX(i).r, fabs(d__1));
L40:
	;
    }
    return ret_val;

/*     End of IZMAX1 */

} /* izmax1_ */

