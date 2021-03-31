#include "Tools.h"

/****************************
Paired-end reads:

R1
--->
_____________________________

                      <------
			   R2

****************************/

int main(int argc, char **argv) {
	
	if( argc != 3)
	{
		std::cerr << "Error usage " << argv[0] << " fileFasta nReads"  << std::endl;
		exit(1);
	}
	
	string fileFasta=argv[1]; 
	dataTypeNSeq n; // numero di Read (#1_F)
	sscanf(argv[2], "%u", &n);

	//Original DA from FASTA file with paired-end reads (R1 and R2) and their RC (R1_RC and R2_RC) in order	
	string fileNameDA=fileFasta+".da";
	FILE * fileDA=fopen(fileNameDA.c_str(), "rb");
	if (fileDA==NULL){
			std::cerr << "Error opening " << fileNameDA << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
	}
	
	//New DA
	string fileNameDA_mod=fileFasta+".da_mod";
	FILE * fileDA_mod=fopen(fileNameDA_mod.c_str(), "wb");
	if (fileDA_mod==NULL){
			std::cerr << "Error opening " << fileNameDA_mod << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
	}
	
	
	//Read original DA
	dataTypeNSeq* DAbuffer=(dataTypeNSeq*)malloc(n*sizeof(dataTypeNSeq));
	int numcharDA=fread(DAbuffer, sizeof(dataTypeNSeq), n, fileDA);
	
	//map DA[i]:
	//2n<=DA[i]<3n (R1_RC) ==> DA[i]=DA[i]-n
	//3n<=DA[i]<4n (R2_RC) ==> DA[i]=DA[i]-3n
	
	while (numcharDA>0) {
		for (dataTypeNChar i=0; i<numcharDA; i++) {
			if (DAbuffer[i]>=2*n && DAbuffer[i]<=(3*n-1)) {
				DAbuffer[i]=DAbuffer[i]-n;
			}
			else if (DAbuffer[i]>=3*n && DAbuffer[i]<=(4*n-1)) {
				DAbuffer[i]=DAbuffer[i]-3*n;
			}	
		}
		
		fwrite(DAbuffer, sizeof(dataTypeNSeq), numcharDA, fileDA_mod);
		// leggi il bloco di n valori successivo
		numcharDA= fread(DAbuffer, sizeof(dataTypeNSeq), n, fileDA);
	}
	
	
	
	fclose(fileDA);
	fclose(fileDA_mod);
	
	//Remove original DA
	if(remove(fileNameDA.c_str()) != 0 )
		cerr << "Error deleting file" << endl;
	//Rename new DA
	if (rename(fileNameDA_mod.c_str(), fileNameDA.c_str()) != 0)
		cerr << "Error renaming file" << endl;
	
	return 0;
}
