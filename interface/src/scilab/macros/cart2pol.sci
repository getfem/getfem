function [r, theta] = cart2pol (x, y)
r = sqrt(x.^2 + y.^2);
old_ieee = ieee();
ieee(2)
theta = atan(y./x);
ieee(old_ieee);
endfunction

