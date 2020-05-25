// Copyright (C) 2009-2020 Yann Colette
// 
//  This file is a part of GetFEM++
// 
//  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
//  under  the  terms  of the  GNU  Lesser General Public License as published
//  by  the  Free Software Foundation;  either version 3 of the License,  or
//  (at your option) any later version along with the GCC Runtime Library
//  Exception either version 3.1 or (at your option) any later version.
//  This program  is  distributed  in  the  hope  that it will be useful,  but
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
//  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
//  License and GCC Runtime Library Exception for more details.
//  You  should  have received a copy of the GNU Lesser General Public License
//  along  with  this program;  if not, write to the Free Software Foundation,
//  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

gf_workspace('clear all');
lines(0);

im  = gf_integ('IM_TRIANGLE(3)');
dim = gf_integ_get(im,'dim');  
gfassert('dim==2');
ise = gf_integ_get(im,'is_exact');
gfassert('~ise');
npt = gf_integ_get(im,'nbpts');
gfassert('npt==[4 2 2 2]');
pts = gf_integ_get(im,'pts');
c   = gf_integ_get(im,'coeffs');
gfassert('size(pts)==[2 10]');
gfassert('size(c)==[1 10]');
im = gf_integ('IM_TRIANGLE(7)');  
c  = gf_integ_get(im,'coeffs');
C  = [0.0267    0.0267    0.0267    0.0386    0.0386    0.0386 ...
0.0386    0.0386    0.0386    0.0878    0.0878    0.0878 ...
-0.0748    0.2460    0.2460    0.4611    0.4611    0.1739 ...
0.1739    0.3261    0.3261    0.1739    0.1739    0.3261    0.3261];
//   C=[0.0386    0.0386    0.0267    0.0267    0.0878    0.0878 ...
//      0.0386    0.0386   -0.0748    0.0878    0.0386    0.0386 ...
//      0.0267    0.2460    0.4611    0.4611    0.2460    0.1739 ...
//      0.3261    0.3261    0.1739    0.1739    0.3261    0.3261    0.1739];
gfassert('norm(c(:)-C(:))<1e-3');
for i=-1:4
  if (i >= 1 & i <= 3) then
    gf_integ_get(im,'face_pts',i);
    gf_integ_get(im,'face_coeffs',i);
  else
    asserterr('gf_integ_get(im,''face_pts'',i)');
    asserterr('gf_integ_get(im,''face_coeffs'',i)');
  end
end
gf_integ_get(im,'char');
asserterr('gf_integ(''IM_TRIANGLE(0)'')');
im  = gf_integ('IM_EXACT_SIMPLEX(3)');
dim = gf_integ_get(im,'dim');  
gfassert('dim==3');
ise = gf_integ_get(im,'is_exact');
gfassert('ise');
asserterr('gf_integ_get(im,''nbpts'')');
gf_integ_get(im,'char');
asserterr('gf_integ(''IM_EXACT_SIMPLEX(0)'')'); // YC: logic exception here: bad parameters ....
