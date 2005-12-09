#!/bin/sh
# http://sources.redhat.com/automake/automake.html#Local-Macros
aclocal -I ./m4
autoheader
autoconf
#pas de ./ devant les noms des makefiles !!!
automake -a --gnu `find . -name Makefile.am | sed -e 's@\./\(.*\)\.am@\1@g'`


