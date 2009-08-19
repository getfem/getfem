function R=axrot_matrix(A, B, theta)
n = (B-A); 
n = n/norm(n);

a=n(1); b=n(2); c=n(3); 

d=sqrt(b^2+c^2);
T=eye(4,4);
T(1:3,4)=-A(:);
Rx=eye(4,4); 

if (norm(n(2:3))>1e-6) then
  Rx(2:3,2:3)=[c -b; b c]/d;
end;

Ry=eye(4,4); 
Ry([1 3],[1 3])=[d -a; a d];

Rz=eye(4,4); 
Rz(1:2,1:2)=[cos(theta) sin(theta); -sin(theta) cos(theta)];

R = inv(T)*inv(Rx)*inv(Ry)*Rz*Ry*Rx*T;
endfunction
