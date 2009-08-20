[X,Y] = meshgrid(-2:0.25:2,-1:0.2:1);
Z = X.* exp(-X.^2 - Y.^2);
[U,V,W] = surfnorm(X,Y,Z);
champ3(X,Y,Z,U,V,W,0.5);
surf(X,Y,Z);
f = gcf();
f.color_map = hsvcolormap(256);


