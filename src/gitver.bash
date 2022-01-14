#!/bin/bash

PATH=$PATH:/usr/bin

git describe --abbrev=6 --dirty --long --always \
  | sed '-es/^\(.*\)$/"\1"/' \
  > gitver.txt

if [ $? != 0 ] ; then
	echo unknown > gitver.txt
fi

cat gitver.txt
