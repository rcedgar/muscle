mkdir -p ../bin

python3 $src/vcxproj_make/vcxproj_make.py --git_hash --openmp --bash --ec2 --binary muscle_ec2 2> make.stderr
rc=$?

echo
echo
echo '=== tail make.stderr ==='
cat make.stderr \
	| grep -v "lto-wrapper: warning: using serial compilation" \
	| grep -v "lto-wrapper: note: see the" \
	| grep -v " warning: Using .dlopen. in statically linked applications" \
	| tail
echo
echo
if [ $rc == 0 ] ; then
	echo SUCCESS
else
	echo ERROR
fi
echo

ls -lh ../bin/muscle_ec2
