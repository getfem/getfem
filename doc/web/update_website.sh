#!/bin/sh

for f in *.php; do
  php $f >  ../../../website/$(basename $f php)html;
done
