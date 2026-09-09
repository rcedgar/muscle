#!/bin/bash

if [ ! -d ../.git ] ; then
  if [ ! -f git_hash.h ] ; then
    echo '#define GIT_HASH "0"' > git_hash.h
  fi
  echo "0" > gitver.txt
  echo "Repo not found, git hash set to zero"
  exit 0
fi

PATH=$PATH:/usr/bin
h=`git describe --abbrev=7 --dirty --long --always`
echo "\"$h\"" > gitver.txt
echo "#define GIT_HASH \"$h\"" > git_hash.h
cat gitver.txt
