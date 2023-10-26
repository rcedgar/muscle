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
  | sed -e 's/\(.*\)/"\1"/' \
  > gitver.txt

cat gitver.txt
