#!/bin/bash

echo -e "Script to install the following tools:\n"
echo -e "- seqtk, toolkit for processing sequences (used to reverse complement sequences);\n"
echo -e "- BCR, tool to build the BWT and related data structures for short read collections (external memory);\n"
echo -e "- eGap, tool to build the BWT and optionally the LCP, DA and SA array (external memory);\n"

mkdir Tools_preprocessing
cd ./Tools_preprocessing

echo -e "\n****Dowloading seqtk...\n"
git clone https://github.com/lh3/seqtk.git;
cd seqtk
echo -e "\n****Compiling seqtk...\n"
make

cd ..
echo -e "\n****Dowloading BCR...\n"
git clone https://github.com/giovannarosone/BCR\_LCP\_GSA
cp ./../BCR_Parameters.h ./BCR\_LCP\_GSA/Parameters.h
cd BCR\_LCP\_GSA
echo -e "\n****Compiling BCR...\n"
make

cd ..
echo -e "\n****Dowloading eGAP...\n"
git clone https://github.com/felipelouza/egap.git
cd egap
echo -e "\n****Compiling eGap...\n"
make


echo -e "\nDone."
