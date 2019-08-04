#!/usr/bin/env bash


javac -Xlint -cp cdk-2.2.jar sgi/*.java 
sleep 5
echo --------------------------------------------------
echo -----------------------start----------------------
echo --------------------------------------------------
echo 5
java -cp cdk-2.2.jar:. -Xms3096m sgi.SGIWithGraphSig 1 5
sleep 30
echo --------------------------------------------------
echo 10
java -cp cdk-2.2.jar:. -Xms3096m sgi.SGIWithGraphSig 1 10
sleep 30
echo --------------------------------------------------
echo 15
java -cp cdk-2.2.jar:. -Xms3096m sgi.SGIWithGraphSig 1 15
sleep 30
echo --------------------------------------------------
echo 20
java -cp cdk-2.2.jar:. -Xms3096m sgi.SGIWithGraphSig 1 20
sleep 30
echo --------------------------------------------------
echo 25
java -cp cdk-2.2.jar:. -Xms3096m sgi.SGIWithGraphSig 1 25
sleep 30
echo -----------------------end------------------------
