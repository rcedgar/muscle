@echo off

if exist c:\cygwin64\bin\bash.exe (
	echo bash found
	c:\cygwin64\bin\bash -c ./gitver.bash
) else (
	echo bash not found
	exit
)
