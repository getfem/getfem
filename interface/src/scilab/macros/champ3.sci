function champ3(x,y,z,fx,fy,fz,c)

X = zeros(2*length(x),1);
Y = zeros(2*length(y),1);
Z = zeros(2*length(z),1);

X(1:2:$) = matrix(x,length(x),1);
X(2:2:$) = matrix(x+fx,length(x),1);
Y(1:2:$) = matrix(y,length(y),1);
Y(2:2:$) = matrix(y+fy,length(y),1);
Z(1:2:$) = matrix(z,length(x),1);
Z(2:2:$) = matrix(z+fz,length(z),1);

xsegs(X,Y,c);
e = gce();
e.arrow_size = 1;
e.data(:,3) = Z;
a = gca();
a.view = '3d';
a.data_bounds(:,3) = [min(e.data(:,3)); max(e.data(:,3))];
endfunction

