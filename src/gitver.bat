@echo off

if exist gitver.txt (
	echo gitver.txt found
) else (
	echo "-" > gitver.txt
)

if exist c:\cygwin64\bin\bash.exe (
	echo bash found
	c:\cygwin64\bin\bash -c ./gitver.bash
) else (
	echo bash not found
	echo 000 > gitver.txt
	exit
)
