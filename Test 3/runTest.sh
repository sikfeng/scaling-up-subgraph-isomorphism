#!/usr/bin/env bash

echo "start compiling"

javac -Xlint -cp cdk-2.2.jar sgi/*.java
sleep 10

dset=(1 2)
gsigs=(1 2 4 8 3 5 6 7 9 10 11)

echo "start running test"
echo "---------------------------------------------------------"
echo ""

echo 1 "direct"
time java -cp cdk-2.2.jar:. -Xms3096m sgi.Direct 1
sleep 60
echo "---------------------------------------------------------"
echo ""

echo 2 "direct"
time java -cp cdk-2.2.jar:. -Xms3096m sgi.Direct 2
sleep 60
echo "---------------------------------------------------------"
echo ""

for gsig in ${gsigs[*]}
do
	for dsetnum in ${dset[*]}
  do
		echo $dsetnum $gsig
		time java -cp cdk-2.2.jar:. -Xms3096m sgi.SGIWithGraphSig $dsetnum $gsig
		sleep 60
		echo "---------------------------------------------------------"
		echo ""
	done
done

