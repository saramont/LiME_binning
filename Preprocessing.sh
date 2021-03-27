#!/bin/bash

########################################################
########################################################
#Script for building DA and LCP from a read FASTA file
########################################################
########################################################

#Please, modify the following commands according to the fasta file paths and names, and the read collection type (paired-end or single-end)

paired=0 #1=paired-end, 0=single-end

PathDataset="./example"
FastaDataset="example.fasta"
FastaDataset_2=""

###############################
Dataset="$(basename "$FastaDataset" .fasta)+RC.fasta"
###############################

###############################
pathseqtk="./Tools_preprocessing/seqtk"
pathBCR="./Tools_preprocessing/BCR_LCP_GSA"
pathEGAP="./Tools_preprocessing/egap"
###############################

#To change the preprocessing tool
BCR=1
EGAP=0

##############################

echo "FastaFile: "$PathDataset/$FastaDataset
#compute the number of reads
nReads=$(grep ">" $PathDataset/$FastaDataset | wc -l)
echo "Number of reads: "$nReads

#Compute reverse complement

if [ $paired -eq 0 ]
then
	FastaDatasetRC="$(basename "$FastaDataset" .fasta)_RC.fasta"
	#reverse complement
	$pathseqtk/seqtk seq -r "$PathDataset/$FastaDataset" > "$PathDataset/$FastaDatasetRC"
	cat $PathDataset/$FastaDataset $PathDataset/$FastaDatasetRC > $PathDataset/$Dataset
	rm 	$PathDataset/$FastaDatasetRC
	#dataset fw+rc
	echo -e "\nDataset: "$PathDataset/$Dataset
else
	FastaDataset_cat="$(basename "$FastaDataset" .fasta)+2.fasta"
    echo -e "\nPaired-end collection\nFastaFile_2: "$PathDataset/$FastaDataset_2
	echo -e "concatenated to: "$PathDataset/$FastaDataset
	#concatenate
    cat $PathDataset/$FastaDataset $PathDataset/$FastaDataset_2 > $PathDataset/$FastaDataset_cat
	FastaDataset_catRC="$(basename "$FastaDataset_cat" .fasta)_RC.fasta"
	#reverse complement
	$pathseqtk/seqtk seq -r "$PathDataset/$FastaDataset_cat" > "$PathDataset/$FastaDataset_catRC"
	cat $PathDataset/$FastaDataset_cat $PathDataset/$FastaDataset_catRC > $PathDataset/$Dataset
	
	rm $PathDataset/$FastaDataset_cat
	rm $PathDataset/$FastaDataset_catRC
	#dataset fw+rc
	echo -e "\nDataset: "$PathDataset/$Dataset   
fi


#####

##Compute DA and LCP

#BCR
if [ $BCR -eq 1 ]
then

    echo "Start BCR..."

    /usr/bin/time -v $pathBCR/BCR_LCP_GSA $PathDataset/$Dataset $Dataset > "BCR_$(basename "$Dataset" .fasta).stdout" 2> "BCR_$(basename "$Dataset" .fasta).stderr"

fi

#EGAP
if [ $EGAP -eq 1 ]
then

    echo "Start eGap..."
	#--mem specify memory assigned to the algorithm in MB. Default is 95% of the available RAM

    /usr/bin/time -v $pathEGAP/eGap $PathDataset/$Dataset --em --mem 4096 --lcp --lbytes 1 --da --dbytes 4 > "eGAP_$(basename "$Dataset" .fasta).stdout" 2> "eGAP_$(basename "$Dataset" .fasta).stderr"
	rm $PathDataset/$Dataset".docs"
	mv $PathDataset/$Dataset".bwt" ./$Dataset".bwt"
	mv $PathDataset/$Dataset".1.lcp" ./$Dataset".lcp"
	mv $PathDataset/$Dataset".4.da" ./$Dataset".da"
fi


if [ $paired -eq 1 ]
then
    ./modifyDA ./DS_dataset/$Dataset $nReads
fi


echo "Done."
####
####
