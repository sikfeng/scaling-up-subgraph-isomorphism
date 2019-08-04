#!/usr/bin/env bash

echo "start compiling"

make -C rdkit-cpp clean
make -C rdkit-cpp direct
make -C rdkit-cpp withgraphsig
javac -Xlint -cp cdk-2.2.jar sgi/*.java
sleep 5

index=(112 991 3505 3599)

echo "start running test"
echo "---------------------------------------------------------"
echo ""

for num in ${index[*]}
do
  echo "rdkit c++ direct" $num
  time ./rdkit-cpp/direct $num
  sleep 60
  echo "---------------------------------------------------------"
  echo ""
  echo "rdkit c++ with graph signatures" $num
  time ./rdkit-cpp/withgraphsig $num
  sleep 60
  echo "---------------------------------------------------------"
  echo ""
  echo "rdkit python direct" $num
  time python direct.py $num
  sleep 60
  echo "---------------------------------------------------------"
  echo ""
  echo "cdk java direct" $num
  time java -cp cdk-2.2.jar:. -Xms3096m sgi.Direct $num
  sleep 60
  echo "---------------------------------------------------------"
  echo ""
  echo "cdk java with graph signatures" $num
  time java -cp cdk-2.2.jar:. -Xms3096m sgi.SGIWithGraphSig $num
  sleep 60
  echo "---------------------------------------------------------"
  echo ""
done
