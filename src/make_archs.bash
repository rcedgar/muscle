mkdir -p ../bin

for arch in x86-64 x86-64-v2 x86-64-v3 x86-64-v4
do
	python3 $src/vcxproj_make/vcxproj_make.py \
		--openmp \
		--bash \
		--arch $arch \
		--binary muscle-$arch \
		2> make.$arch.stderr

	rc=$?

	echo
	echo
	echo '=== tail make.$arch.stderr ==='
	cat make.$arch.stderr \
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

	ls -lh ../bin/muscle-$arch
done
