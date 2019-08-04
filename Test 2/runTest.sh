#!/usr/bin/env bash

echo "start compiling"

javac -Xlint -cp cdk-2.2.jar benchmark.java
sleep 10

echo cdk 1
time java -cp cdk-2.2.jar:. -Xms3096m benchmark 1
sleep 60

echo igraph 1
time python benchmark.py 1
sleep 60

echo cdk 2
time java -cp cdk-2.2.jar:. -Xms3096m benchmark 2
sleep 60

echo python 2
time python benchmark.py 2
