# LiME_binning

LiME_binning is a novel lightweight alignment-free and assembly-free framework for metagenomic binning that is combinatorial by nature and allows us to use little internal memory. 

Let *S* be a large collection of biological sequences comprising reads.

It takes in input:
- the longest common prefix array (lcp) of collection *S*;
- the document array (da) of collection *S*.


### Install

```sh
git clone https://github.com/saramont/LiME_binning
cd LiME_binning
```

### Compile

```sh
make
```

### Run
```sh
./
```

### Examples
```sh
./
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
