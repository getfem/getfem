function varargout=gf_colormap(name),
// function  c=gf_colormap(name)
//   return a colormap, or change the current colormap.
//   name can be: 'tripod', 'chouette', 'froid', 'tank'
//   or 'earth'.

if (name=='tripod') then
  r=[0.7 .7 .7]; 
  l = r($,:); 
  s=63; s1=20; s2=25; s3=48;s4=55; 
  for i=1:s, 
    c1 = max(min((i-s1)/(s2-s1),1),0);
    c2 = max(min((i-s3)/(s4-s3),1),0); 
    r($+1,:)=(1-c2)*((1-c1)*l + c1*[1 0 0]) + c2*[1 .8 .2]; 
  end;
elseif (name=='chouette') then
    gg = [.8 1.0 .8; 
          .7  .9 .4;
          .3  .8 .2;
          .1  .7 .4;
          .2 0.7 1.0;
          .3 0.3 1.0;
         1.0  .8  .1;
         1.0  .6  .1;
         1.0  .45 .1;
         1.0 0.3  .1];
    r = matrix(repmat(gg',6,1),3,60)'; 
elseif (name=='froid') then
    gg = [ .8 1.0  .8; 
           .7  .9  .4;
           .3  .8  .2;
           .1  .7  .4;
           .2 0.7 1.0;
           .3 0.3 1.0];
    r = matrix(repmat(gg',10,1),3,60)';     
elseif (name=='tank') then
  r=[0  0 1; 0    .5 1; 0 1    .5; 
     0  1 0;  .5 1   0; 1  .5 0; 
     1 .4 0; 1 0     0; 1  .2 0; 
     1 .4 0; 1    .6 0; 1  .8 0];
    r = matrix(repmat(r',5,1),3,60)';     
elseif (name=='earth') then
  r=[252 233  79; //	Butter 1
     247 222  30;
     237 212   0; //	Butter 2
     216 180   0;
     196 160   0; //	Butter 3
     138 226  52; //	Chameleon 1
     115 210  22; //	Chameleon 2
     78 154   6];
  r = matrix(r'/255, 8,1);
  r = matrix(r,3,length(r)/3)';     
else
  error('wrong colormap');
end
if (nargout) then
  varargout=list(r); // YC: to be verified
else
  f = gcf();
  f.color_map = r;
end

c=[252 233  79; //	Butter 1
   237 212   0; //	Butter 2
   196 160   0; //	Butter 3
   138 226  52; //	Chameleon 1
   115 210  22; //	Chameleon 2
   78 154   6; //	Chameleon 3
   252 175  62; //	Orange 1
   245 121   0; //	Orange 2
   206  92   0; //	Orange 3
   114 159 207; //	Sky Blue 1
   114 159 207; //	Sky Blue 1
   52 101 164; //	Sky Blue 2
   52 101 164; //	Sky Blue 2
   32  74 135; //	Sky Blue 3
   32  74 135; //	Sky Blue 3
   173 127 168; //	Plum 1
   173 127 168; //	Plum 1
   173 127 168; //	Plum 1
   117  80 123; //	Plum 2
   117  80 123; //	Plum 2
   117  80 123; //	Plum 2
   92  53 102; //	Plum 3
   92  53 102; //	Plum 3
   92  53 102; //	Plum 3
   233 185 110; //	Chocolate 1
   233 185 110; //	Chocolate 1
   233 185 110; //	Chocolate 1
   233 185 110; //	Chocolate 1
   193 125  17; //	Chocolate 2
   193 125  17; //	Chocolate 2
   193 125  17; //	Chocolate 2
   193 125  17; //	Chocolate 2
   143  89   2; //	Chocolate 3
   143  89   2; //	Chocolate 3
   143  89   2; //	Chocolate 3
   143  89   2; //	Chocolate 3
   239  41  41; //	Scarlet Red 1
   239  41  41; //	Scarlet Red 1
   239  41  41; //	Scarlet Red 1
   239  41  41; //	Scarlet Red 1
   239  41  41; //	Scarlet Red 1
   204   0   0; //	Scarlet Red 2
   204   0   0; //	Scarlet Red 2
   204   0   0; //	Scarlet Red 2
   204   0   0; //	Scarlet Red 2
   204   0   0; //	Scarlet Red 2
   164   0   0; //	Scarlet Red 3
   164   0   0; //	Scarlet Red 3
   164   0   0; //	Scarlet Red 3
   164   0   0; //	Scarlet Red 3
   164   0   0]; //	Scarlet Red 3
//     238 238 236]; //	Untitled

c=c/255;
endfunction

