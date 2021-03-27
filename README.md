# LiME_binning

LiME_binning is a novel lightweight alignment-free and assembly-free framework for metagenomic binning that is combinatorial by nature and allows us to use little internal memory. 

Let *S* be a large collection of biological sequences comprising reads.

It takes in input:
- the longest common prefix array (lcp) of collection *S*;
- the document array (da) of collection *S*.


The underlying method can be summarized in two main steps:

1. By reading lcp(S), we detect alpha-clusters, i.e., blocks of ebwt(S) symbols whose associated suffixes share a common context of a minimum length alpha;

2. By reading da(S), we analyze alpha-clusters to evaluate a degree of similarity for any pair of two reads in S. 
Then, we group reads together according to the similarity degree and a minimum threshold value t. 



### Install

```sh
git clone https://github.com/saramont/LiME_binning
cd LiME_binning
make
```

###Preprocessing step

There are several tools available for computing the required data structures (lcp, da) for S.

For short reads, one could use BCR [https://github.com/giovannarosone/BCR_LCP_GSA], or eGAP [https://github.com/felipelouza/egap], for instance.
To install BCR and eGap for the preprocessing, one could run
 ```sh
./Install_preprocessing_tools.sh
```
 
To build up the data structures for the collection S, one could run
 ```sh
./Preprocessing.sh
```
by first setting in it the fasta file name and its path.

In addition, if S is a paired-end collection, one could set paired=1 in the script Preprocessing.sh, before running it.

### Run
The steps of the method are accomplished by running:

(1) ClusterLCP with input parameters: name of the fasta file, number of reads in S, alpha;
```sh
./ClusterLCP fileName numReads alpha
```
(2) BinningDA with input parameters: name of the fasta file, minimum threshold similarity value (t), number of matrix rows in main memory (value in [1,numReads]);
```sh
./BinningDA fileName t max_rows
```
The output file of BinningDA is a text file of numReads lines where the binning results are stored according to the format:
the i-th line contains the ID-number group which the i-th read is associated with.

One may use the script LiME_binning.sh to run both steps sequentially. It takes as input:
- fileName, the fasta file name;
- numReads, the number of reads in S;
- alpha, minimum LCP value in a cluster;
- t, the minimum threshold similarity value.

One may reduce the internal memory usage by setting the parameter max_rows (max_rows is value in [1,numReads]).
Indeed, such parameter limits the number of matrix rows in main memory to max_rows. (Default: numReads).

```sh
./LiMEbinning fileName numReads alpha t
```


### Examples
After running Install_tools_preprocessing.sh and Preprocessing.sh (on input example.fasta):
```sh
./LiMEbinning example+RC.fasta 12 16 20
```



## References:

Degree thesis
```sh
@thesis{thesisMontemaggi,
  author       = {Sara Montemaggi}, 
  title        = {Nuovo approccio per il binning di frammenti di DNA: teoria ed esperimenti},
  school       = {University of Pisa},
  note         = {Supervisors: Giovanna Rosone and Veronica Guerrini}
}
```

---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html
