#!/bin/bash -e

PATH=$PATH:/usr/bin

git describe --abbrev=6 --dirty --long --always \
  | sed '-es/^\(.*\)$/"\1"/' \
  > gitver.txt

cat gitver.txt
