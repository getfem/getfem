#include "slu_Cnames.h"

int lsame_(char *ca, char *cb)
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

    LSAME returns .TRUE. if CA is the same letter as CB regardless of case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
*/  

    /* System generated locals */
    int ret_val;
    
    /* Local variables */
    int inta, intb, zcode;

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

    /* Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

    /* Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {
	/* ASCII is assumed - ZCODE is the ASCII code of either lower or   
          upper case 'Z'. */
	if (inta >= 97 && inta <= 122) inta += -32;
	if (intb >= 97 && intb <= 122) intb += -32;

    } else if (zcode == 233 || zcode == 169) {
	/* EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or   
          upper case 'Z'. */
	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169)
	    inta += 64;
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169)
	    intb += 64;
    } else if (zcode == 218 || zcode == 250) {
	/* ASCII is assumed, on Prime machines - ZCODE is the ASCII code   
          plus 128 of either lower or upper case 'Z'. */
	if (inta >= 225 && inta <= 250) inta += -32;
	if (intb >= 225 && intb <= 250) intb += -32;
    }
    ret_val = inta == intb;
    return ret_val;
    
} /* lsame_ */
