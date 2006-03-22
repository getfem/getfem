#!/bin/sh

# update macros for linking to doxygen doc.
rm filelinks.tex;
cd ../../src 
for i in *.h *.cc; do
    echo "doing ${i}";
    k=$(echo $i | sed -e 's/[._]//g' -e 's/2/two/g' -e 's/1/one/g' -e 's/0/zero/g');
    j=$(echo $i | sed -e 's/_/__/g' -e 's/\./_8/' );
    i=$(echo $i | sed -e 's/_/\\_/g');
    echo "\newcommand{\\${k}}{\\doxfilename{$i}{$j}\xspace}" >> ../doc/userdoc/filelinks.tex;
done;
