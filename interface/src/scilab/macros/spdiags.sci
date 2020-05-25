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

function [A_out,d_out] = spdiags(B_in,d_in,n_row,n_col)

// B = spdiags(A) extracts all nonzero diagonals from the m-by-n matrix A. B is a min(m,n)-by-p matrix whose columns are the p nonzero diagonals of A.
// [B,d] = spdiags(A) returns a vector d of length p, whose integer components specify the diagonals in A.
// B = spdiags(A,d) extracts the diagonals specified by d.
// A = spdiags(B,d,A) replaces the diagonals specified by d with the columns of B. The output is sparse.
// OK - A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.

// spdiags(matrix(1:12, 4, 3), [-1 0 1], 5, 4)
// result: 5 10  0  0
//         1  6 11  0
//         0  2  7 12
//         0  0  3  8
//         0  0  0  4

[nargout,nargin] = argn();

A_out = sparse([]);
d_out = [];
diagonal_extract = %F;

// Diagonal extraction
if (nargin==1) then
  n_row = size(B_in,1);
  n_col = size(B_in,2);
  n_min = min(n_row,n_col);
  d_in  = -(n_min-1):(n_min-1); 
  A_out = spzeros(n_min,length(d_in));
  diagonal_extract = %T;
end

if (nargin==2) then
  n_row = size(B_in,1);
  n_col = size(B_in,2);
  n_min = min(n_row,n_col);
  A_out = spzeros(n_min,length(d_in));
  diagonal_extract = %T;
end

// Diagonal matrix creation
if (nargin==3) then
  A_out = sparse(n_row);
  n_row = size(A_out,1);
  n_col = size(A_out,2);
  diagonal_extract = %F;
end

if (nargin==4) then
  A_out = spzeros(n_row,n_col);
  diagonal_extract = %F;
end

if (nargin>4) | (nargin<1) then
  error('1 to 4 argument required');
end

if (max(d_in)>=n_col) | (min(d_in)<=-n_row) then
  error('diagonal index not in range');
end

n_min = min(n_row,n_col);

if (diagonal_extract) then
  for i=1:length(d_in)
    printf('d(%d) = %d\n', i, d_in(i));
    disp(size(A_out))
    if (d_in(i)>0) then
      A_out(1:n_min-d_in(i),i) = B_in(sub2ind(size(B_in),1:n_min-d_in(i),1+d_in(i):n_min));
    elseif (d_in(i)<0) then
      A_out(1-d_in(i):n_min,i) = B_in(sub2ind(size(B_in),1-d_in(i):n_min,1:n_min+d_in(i)));
    else
      A_out(1:n_min,i) = B_in(sub2ind(size(B_in),1:n_min,1:n_min));
    end
  end
  if (nargout==2) then
    d_out = [];
    for i=size(A_out,2):-1:1
      if and(A_out(:,i)==0) then
        A_out(:,i) = [];
      else
        d_out = [i d_out];
      end
    end
    d_out = d_out - n_min;
  end
  
else
  for i=1:length(d_in)
    if (d_in(i)>0) then
      n_row_start = 1;
      n_row_end   = n_min;
      n_col_start = d_in(i)+1;
      n_col_end   = min(n_min+d_in(i)+1, n_col);
      mat_size    = [n_row_end - n_row_start + 1 n_col_end - n_col_start + 1]
      b_start     = 1;
      b_end       = min(size(B_in,1),min(mat_size));
    elseif (d_in(i)<0) then
      n_row_start = -d_in(i)+1;
      n_row_end   = min(n_min-d_in(i)+1, n_row);
      n_col_start = 1;
      n_col_end   = n_min;
      mat_size    = [n_row_end - n_row_start + 1 n_col_end - n_col_start + 1]
      b_start     = -d_in(i)+1;
      b_end       = size(B_in,1);
    else
      n_row_start = 1;
      n_col_start = 1;
      n_row_end   = n_min;
      n_col_end   = n_min;
      mat_size    = [n_row_end - n_row_start + 1 n_col_end - n_col_start + 1]
      b_start     = 1;
      b_end       = n_min;
    end

    A_out(n_row_start:n_row_end,n_col_start:n_col_end) = A_out(n_row_start:n_row_end,n_col_start:n_col_end) + ...
                                    sparse([1:b_end-b_start+1;1:b_end-b_start+1]', ...
                                           B_in(b_start:b_end,i), ...
                                          [n_row_end - n_row_start + 1 n_col_end - n_col_start + 1]);
  end
end
endfunction
