function [r, theta] = cart2pol (x, y)
r = sqrt(x.^2 + y.^2);
theta = atan (y./x);
endfunction

