#include "slu_Cnames.h"
#include "slu_dcomplex.h"

double dzsum1_(int *n, doublecomplex *cx, int *incx)
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

    DZSUM1 takes the sum of the absolute values of a complex   
    vector and returns a double precision result.   

    Based on DZASUM from the Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with ZLACON.   

    Arguments   
    =========   

    N       (input) INT   
            The number of elements in the vector CX.   

    CX      (input) COMPLEX*16 array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INT   
            The spacing between successive values of CX.  INCX > 0.   

    ===================================================================== 
*/  

    /* Builtin functions */
    double z_abs(doublecomplex *);
    
    /* Local variables */
    int i, nincx;
    double stemp;


#define CX(I) cx[(I)-1]

    stemp = 0.;
    if (*n <= 0) {
	return stemp;
    }
    if (*incx == 1) {
	goto L20;
    }

    /*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {

	/*        NEXT LINE MODIFIED. */

	stemp += z_abs(&CX(i));
/* L10: */
    }
    
    return stemp;

    /*     CODE FOR INCREMENT EQUAL TO 1 */

L20:
    for (i = 1; i <= *n; ++i) {

	/*        NEXT LINE MODIFIED. */

	stemp += z_abs(&CX(i));
/* L30: */
    }
    
    return stemp;

    /*     End of DZSUM1 */

} /* dzsum1_ */

