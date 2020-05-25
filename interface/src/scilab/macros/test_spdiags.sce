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

// Test 1

A = [0 5 0 10  0  0;...
     0 0 6  0 11  0;...
     3 0 0  7  0 12;...
     1 4 0  0  8  0;...
     0 2 5  0  0  9];
   
[B, d] = spdiags(A);

B_ref = [0 0 5 10; ...
         0 0 6 11; ...
         0 3 7 12; ...
         1 4 8  0; ...
         2 5 9  0];
d_ref = [-3 -2 1 3];

if (~and(B==B_ref) | ~and(d==dref)) then
  printf('error in test1\n');
end

// Test 2

n = 10;
e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
A = spdiags(abs(-(n-1)/2:(n-1)/2)',0,A);
B = spdiags(A);

// Test 3

A = [11    0   13    0
      0   22    0   24
      0    0   33    0
     41    0    0   44
      0   52    0    0
      0    0   63    0
      0    0    0   74];
      
[B,d] = spdiags(A);

d_ref = [-3 0 2];

B_ref = [41 11  0; ...
         52 22  0; ...
         63 33 13; ...
         74 44 24]

if (~and(B==B_ref) | ~and(d==dref)) then
  printf('error in test 3\n');
end

// Test 4

B = [1:6]' .*. ones(1,7);

//B =
//
//    1  1  1  1  1  1  1
//    2  2  2  2  2  2  2
//    3  3  3  3  3  3  3
//    4  4  4  4  4  4  4
//    5  5  5  5  5  5  5
//    6  6  6  6  6  6  6

d = [-4 -2 -1 0 3 4 5];
A = spdiags(B,d,6,6);

A_ref = [1 0 0 4 5 6; ...
         1 2 0 0 5 6; ...
         1 2 3 0 0 6; ...
         0 2 3 4 0 0; ...
         1 0 3 4 5 0; ...
         0 2 0 4 5 6];

if (~and(A==A_ref)) then
  printf('error in test 4\n');
end

