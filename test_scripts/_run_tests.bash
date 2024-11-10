#!/bin/bash -e

cd ../src
rm -rf o/ ../bin/muscle*
chmod +x ./build_linux.bash
./build_linux.bash
cd ../test_scripts

rm -rf ../test_output ../test_logs
mkdir -p ../test_output ../test_logs ../test_results

log=../test_output/TEST_LOG.txt

date=`date "+%Y-%m-%d/%H:%M:%S"`
ver=`../bin/muscle --version | tr -d ' \n\r' | sed "-es/Built.*//"`
echo $date $ver STARTED >> $log
git status >> $log

echo STARTED `date` >> $log

./run_BB_seqs.bash
./run_BB_structs.bash
./run_rdrp_seqs.bash
./run_rdrp_structs.bash

./qscore_BB_seqs.bash
./qscore_BB_structs.bash
./qscore_rdrp.bash

python3 ./check_logs.py
python3 ./check_results.py
python3 ./update_success_list.py $ver $date

echo COMPLETED $date >> $log

echo $date $ver SUCCESS >> $log
