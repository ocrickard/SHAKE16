#!/bin/sh

#  test.sh
#  SHAKE16
#
#  Created by Oliver Rickard on 7/9/16.
#  Copyright © 2016 Oliver Rickard. All rights reserved.

rm -rf output/test
mkdir -p output/test

cd test-data

SHAKE16 INP.DAT ../output/test/output1.txt ../output/test/output2.txt

cd ..

file1diff=$(diff -u test-data/output1.txt output/test/output1.txt)
file2diff=$(diff -u test-data/output2.txt output/test/output2.txt)

if [ "$file1diff" == "" ]
then
  echo "output file 1 identical, test passed"
else
  echo "output file 1 did not match, test failed. diff:"
  echo $file1diff
fi

if [ "$file2diff" == "" ]
then
echo "output file 2 identical, test passed"
else
echo "output file 2 did not match, test failed. diff:"
echo $file1diff
fi
