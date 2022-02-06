#!/bin/bash

if [ ! -d ../.git ] ; then
  if [ ! -f gitver.txt ] ; then
    echo "0" > gitver.txt
  fi
  echo "Repo not found, git hash set to zero"
  exit 0
fi

PATH=$PATH:/usr/bin

git describe --abbrev=6 --dirty --long --always \
  > gitver.tmp

sed -i '-es/"//g' gitver.tmp

echo \"`cat gitver.tmp`\" > gitver.txt

rm -f gitver.tmp

cat gitver.txt
