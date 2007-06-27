#!/bin/sh
IN=$1
OUT=$2
bf=$(basename $IN)
TMP=$(mktemp /tmp/${bf}.XXXXXX) || exit 1
#TMP=./demo_laplacian.tmp
echo '\\begin{mcode}' > $TMP
sed -e 's/{/\\{/g' -e 's/}/\\}/g' < $IN >> $TMP
echo '\\end{mcode}' >> $TMP

perl $(dirname $0)/latexize_mcode.pl < $TMP > $OUT || ( rm $OUT; exit 1)
