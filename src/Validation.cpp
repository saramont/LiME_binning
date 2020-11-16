#include "Tools.h"


//arguments: file_name num_genomes num_read_g1 num_read_g2 ...

int checkGenome(dataTypeNSeq readID, int numGenomes, dataTypeNSeq *readsPerGenome) {
	dataTypeNSeq start=0;
	dataTypeNSeq end=readsPerGenome[0]-1;
	
	for (int i=0; i<numGenomes; i++) {
		if (readID>=start && readID<= end) {
			return i;
		}
		start=end+1;
		end=start+readsPerGenome[i+1]-1;
	}
	return -1;
}


int main(int argc, char **argv) {
	
	//Read arguments
	if (argc < 4) {	
		std::cerr << "Error usage: must have at least 3 arguments" << endl; 
		exit(1);
	}
	
	string fnOutSets=argv[1];
	
	int numGenomes;
	sscanf(argv[2], "%d", &numGenomes);
	if (argc != numGenomes+3) {
		std::cerr << "Error usage: incorrect number of arguments" << endl; 
		exit(1);
	}
	
	dataTypeNSeq *readsPerGenome = (dataTypeNSeq *)malloc(numGenomes*sizeof(dataTypeNSeq));
	dataTypeNSeq numReads=0;
	for (int i=0; i<numGenomes; i++) {
		sscanf(argv[3+i], "%d", &readsPerGenome[i]); 
		numReads=numReads+readsPerGenome[i];
	}
	cout << "Total number of reads: " << numReads << endl;
	
	float thresh=numReads*0.1/100;
	
	cout << "Thresh: " << thresh << endl;
	
	//Read output file with set assignments
	FILE * fpOutSets;
	fpOutSets=fopen(fnOutSets.c_str(), "r");
	if(fpOutSets==NULL) {
		std::cerr << "Error opening " << fnOutSets << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
	//readToSet[i] = setID of the i-th read
	dataTypeNSeq *readToSet = (dataTypeNSeq *)malloc(numReads*sizeof(dataTypeNSeq));
	for (dataTypeNSeq i=0; i<numReads; i++) {
		dataTypeNSeq readID;
		//fscanf(fpOutSets, "%u", &readID);
		fscanf(fpOutSets, "%u", &readToSet[i]);
	}
	fclose(fpOutSets);
	
	//Get sets dimensions to discard sets with less than 100 reads
	map<dataTypeNSeq, dataTypeNSeq> setsDim; //key = setID, value = number of elements in set setID
	for (dataTypeNSeq i=0; i<numReads; i++) {
		setsDim[readToSet[i]]++;	
	}
	
	//Calculate new total number of sets and updated sets IDs
	map<dataTypeNSeq, dataTypeNSeq> oldToNewID;
	dataTypeNSeq nextID=0;
	dataTypeNSeq totNonGrouped=0;
	dataTypeNSeq numSets=0;
	for (map<dataTypeNSeq, dataTypeNSeq>::iterator it=setsDim.begin(); it!=setsDim.end(); ++it) {
		if (it->second >= thresh) {
			numSets++;
			oldToNewID.insert(pair<dataTypeNSeq,dataTypeNSeq>(it->first, nextID));
			nextID++;
		}
		else 
			totNonGrouped=totNonGrouped+it->second;
	}
	
	cout << "Number of sets: " << numSets << endl;
	
	//Construct matrix
	dataTypeNSeq **A=(dataTypeNSeq **)malloc(numSets*sizeof(dataTypeNSeq *));
	for (dataTypeNSeq i=0; i<numSets; i++) {
		A[i] = (dataTypeNSeq *)malloc(numGenomes*sizeof(dataTypeNSeq));
		for (int j=0; j<numGenomes; j++) {
			A[i][j]=0;
		}
	}
	for (dataTypeNSeq i=0; i<numReads; i++) {
		if (setsDim[readToSet[i]] >= thresh) {  
			//update corresponding value
			int genome=checkGenome(i, numGenomes, readsPerGenome);
			A[oldToNewID[readToSet[i]]][genome]++;
		}
	}	
	
	//Calculate precision, recall and F-measure
	dataTypeNSeq sumRowMaxes=0;
	dataTypeNSeq sumColMaxes=0;
	dataTypeNSeq sumAij=numReads-totNonGrouped;
	dataTypeNSeq *rowMaxes = (dataTypeNSeq *)malloc(numSets*sizeof(dataTypeNSeq));
	for (dataTypeNSeq i=0; i<numSets; i++) {
		rowMaxes[i]=0;
	}
	dataTypeNSeq *colMaxes = (dataTypeNSeq *)malloc(numGenomes*sizeof(dataTypeNSeq));
	for (int i=0; i<numGenomes; i++) {
		colMaxes[i]=0;
	}
	for (dataTypeNSeq i=0; i<numSets; i++) {
		for (int j=0; j<numGenomes; j++) {
			if (A[i][j] > rowMaxes[i]) rowMaxes[i]=A[i][j];
			if (A[i][j] > colMaxes[j]) colMaxes[j]=A[i][j];
		}
	}
	
	for (dataTypeNSeq i=0; i<numSets; i++) {
		sumRowMaxes=sumRowMaxes+rowMaxes[i];
	}
	for (int i=0; i<numGenomes; i++) {
		sumColMaxes=sumColMaxes+colMaxes[i];
	}
	
	double precision=(double)sumRowMaxes/(double)sumAij;
	double recall=(double)sumColMaxes/(double)(sumAij+totNonGrouped);
	double F_measure=(2*precision*recall)/(precision+recall);
	
	cout << "Precision: " << precision << endl;
	cout << "Recall: " << recall << endl;
	cout << "F-measure: " << F_measure << endl;
	

	return 0;
}
