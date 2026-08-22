#!/bin/bash -e

if [ ! -s ./_run_tests.bash ] ; then
	echo "ERROR -- must run in test_scripts/ directory"
	exit 1
fi

./_run_tests.bash

if [ $? != 0 ] ; then
	echo === FAILED ===
fi

tail ../test_output/TEST_LOG.txt
