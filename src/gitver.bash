#!/bin/bash

PATH=$PATH:/usr/bin

git describe --abbrev=6 --dirty --long --always \
  > gitver.tmp

sed -i '-es/"//g' gitver.tmp

echo \"`cat gitver.tmp`\" > gitver.txt

rm -f gitver.tmp

cat gitver.txt
