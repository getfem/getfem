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

if (nargin==4) then
  A_out = spzeros(n_row,n_col);
  diagonal_extract = %F;
end

if (nargin>4) | (nargin<1) then
  error('1 to 4 argument required');
end

//if (max(d_in)>=n_col) | (min(d_in)<=-n_row) then
//  error('diagonal index not in range');
//end

horiz_matr = (size(B_in,1)<size(B_in,2));
n_min = min(n_row,n_col);

if (diagonal_extract) then
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
