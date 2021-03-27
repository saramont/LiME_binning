#!/bin/bash

#Input
FastaFilename=$1
numReads=$2

#Parameters
alpha=$3
t=$4

#Set max_rows in [1,numReads) to limit memory usage 
#note that the process may slow down heavily for very low values of max_rows
#binningDA requires max_rows*2*numReads Bytes in memory
max_rows=$numReads


echo "LiMEbinning:"
echo -e "\talpha = "$alpha
echo -e "\tminimum similarity value = "$t
echo -e "\tnumber of matrix rows in main memory = "$max_rows


#Choose the steps you want to execute
step1=1   #=1 to execute ClusterLCP.
step2=1   #=1 to execute BinningDA.


#First step
if [ $step1 -eq 1 ]
then
	echo "Running ClusterLCP..."
	/usr/bin/time -v ./ClusterLCP $FastaFilename $numReads $alpha > "ClusterLCP_"$FastaFilename".stdout" 2> "ClusterLCP_"$FastaFilename".stderr"

fi

#Second step
if [ $step2 -eq 1 ]
then
    echo "Running BinningDA..."
	/usr/bin/time -v ./BinningDA $FastaFilename $t $max_rows > "BinningDA_"$FastaFilename".stdout" 2> "BinningDA_"$FastaFilename".stderr"

fi
