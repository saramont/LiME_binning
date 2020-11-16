#include "Tools.h"

/*
Step 2. analyzes alpha-clusters in order to evaluate a degree of similarity between any readsin the collection.

As preprocessing, it need to have fileFasta.da computed.

Input: fileFasta, minimum similarity value t, maximum number of rows in memory, threads
 
Output: Resume on sets created and file with sets associations. Optionally binary files with all non-zeros similarity values
 
*/


#ifndef SMALL
	#define SMALL 0
#endif

#ifndef OMP
	#define OMP 1
#endif

#ifndef STEP2
	#define STEP2 0
#endif

#ifndef SAVE_MAT
	#define SAVE_MAT 0
#endif

#ifndef SQUARE
	#define SQUARE 1
#endif


#if STEP2
void Analysis_and_updating(std::vector<dataTypeNSeq>::iterator &lowBounds, std::vector<dataTypeNSeq>::iterator &upBounds, vector<dataTypeSetSim> &read_symb, dataTypeNSeq n, dataTypeSimsMatS2 &simVals, dataTypeNSeq firstRow, dataTypeNSeq lastRow)
#else
void Analysis_and_updating(std::vector<dataTypeNSeq>::iterator &lowBounds, std::vector<dataTypeNSeq>::iterator &upBounds, vector<dataTypeSim> &read_symb, dataTypeNSeq n, dataTypeSimsMat &simVals, dataTypeNSeq firstRow, dataTypeNSeq lastRow)
#endif
{	
	
    vector<dataTypeNSeq> id_reads;
    
    dataTypeNSeq da, da_next, idxStart, idxEnd;
    bool foundStart = false;
    bool foundEnd = false;
	//Reads
    while ( lowBounds < upBounds ){ 
		// read_symb[i] = quanti simboli, nell'alpha-cluster che sto al momento analizzando, vengono dal read i
		da = *lowBounds;
        da_next = *lowBounds;
        read_symb[da]=0; 
 
        id_reads.push_back(da); // id_reads contiene gli indici ordinati!
        if (!foundStart && da>=firstRow && da<=lastRow) { 
			foundStart = true;
			foundEnd = true;
			idxStart = id_reads.size()-1;
			idxEnd = id_reads.size()-1;
		}
		else if (foundStart && da>=firstRow && da<=lastRow) {
			idxEnd = id_reads.size()-1;
		}
        while ( ( da == da_next ) && ( lowBounds < upBounds ) ) 
        {
            read_symb[da]++;
            lowBounds++;
            da = da_next;
			da_next = *lowBounds;
        }
    }
    // a fine while controllare che la dimensione di id_reads sia maggiore di 1, se è 1 questo
    // cluster non mi interessa e non devo fare nulla --> return
    // controlla inoltre che in questo cluster sia presente l'indice della riga che sto analizzando,
    // altrimenti è inutile
    if (id_reads.size() == 1 || !foundStart) return;

    #if STEP2
		dataTypeSetSim t=0;
    #else
		dataTypeSim t=0;
	#endif
	// se arrivo qua sono sicura di avere più di un read e di avere almeno un read delle righe che sto analizzando
	for (dataTypeNSeq i=idxStart; i<=idxEnd; i++) //range over reads 
    {
		#if SQUARE
		for (dataTypeNSeq j=0; j<id_reads.size(); j++)
		#else
		for (dataTypeNSeq j=i+1; j<id_reads.size(); j++)
        #endif
		{
			#if SQUARE
			if ((id_reads[j]%n) != id_reads[i]) // non salvo le similarità read_i con read_i (o con complementare read_i)
			#else
			if ((id_reads[j]%n) > id_reads[i]) 
			#endif
			{
            //Take the minimum
				if (read_symb[id_reads[i]] > read_symb[id_reads[j]] ) 
					t= read_symb[id_reads[j]]; 
				else
					t= read_symb[id_reads[i]]; 

				//Update similarity in current row
				if(t>0){
					#if OMP
					#pragma omp critical
					#endif
					simVals[id_reads[i]-firstRow][id_reads[j]] += t;
				}
			}//end-if	
		}//end-for
    }//end-for
	return;
}

#if STEP2
dataTypeNChar clusterAnalyze(dataTypeNChar chk, FILE *InFileCluster, FILE *InDA, dataTypeNSeq n, dataTypeSimsMatS2 &simVals, dataTypeNSeq firstRow, dataTypeNSeq lastRow)
#else
dataTypeNChar clusterAnalyze(dataTypeNChar chk, FILE *InFileCluster, FILE *InDA, dataTypeNSeq n, dataTypeSimsMat &simVals, dataTypeNSeq firstRow, dataTypeNSeq lastRow) // tolto nRef dagli argomenti
#endif
{
	dataTypeNChar counter=0;
    dataTypeNChar numcharCluster, numcharDA;
	ElementCluster* clusterbuffer= new ElementCluster[BUFFERCLUSTER];
    dataTypeNSeq* DAbuffer= new dataTypeNSeq[sizeMaxBuf]; //65536
	
	
	vector<dataTypeNSeq> elebuffer;
	std::vector<dataTypeNSeq>::iterator lowBounds, upBounds; // pos_nReadInCluster non serve più
	
	elebuffer.resize(sizeMaxBuf);
    
	#if STEP2
		vector<dataTypeSetSim> read_symb;
		read_symb.resize(n);
	#else
		vector<dataTypeSim> read_symb;
		read_symb.resize(2*n);
	#endif
    
    //Read InFileCluster
    if(BUFFERCLUSTER<chk)
		numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),BUFFERCLUSTER,InFileCluster);
    else
        numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),chk,InFileCluster);

        
    while(numcharCluster>0)
    {
        for (dataTypeNChar c=0;c<numcharCluster; c++) // scorre su tutte le coppie, che rappresentano i cluster, che sono state lette
		{
            //Read cluster DA and EBWT
			fseek(InDA, (clusterbuffer[c].pStart)*sizeof(dataTypeNSeq), SEEK_SET);
			numcharDA=fread(&elebuffer[0],sizeof(dataTypeNSeq),clusterbuffer[c].len, InDA);
			
			//Sort
			std::sort (elebuffer.begin(), elebuffer.begin()+ numcharDA);//w.r.t. da
			                
			//Set lowBounds and upBounds
			lowBounds=elebuffer.begin();
			upBounds=elebuffer.begin()+ numcharDA;
            Analysis_and_updating(lowBounds, upBounds, read_symb, n, simVals, firstRow, lastRow); 
            
        }//end-for clusters
        counter+=numcharCluster;
        chk-=numcharCluster;
        if(chk>0 && BUFFERCLUSTER<chk)
            numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),BUFFERCLUSTER,InFileCluster);
        else
            numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),chk,InFileCluster);
    }//end-while
    elebuffer.shrink_to_fit();
	read_symb.shrink_to_fit();
	
	return counter;
}


#if STEP2
void setAnalysis(dataTypeSimsMatS2 &simVals, dataTypeNSeq *reads_to_set, dataTypeNSeq nRows, dataTypeNSeq firstRow, float t, dataTypeNSeq n, dataTypeNSeq &next_set_id, dataTypeNSeq &tot_set, dataTypeNSeq *sizes) 
#else
void setAnalysis(dataTypeSimsMat &simVals, dataTypeNSeq *reads_to_set, dataTypeNSeq nRows, dataTypeNSeq firstRow, dataTypeSim t, dataTypeNSeq n, dataTypeNSeq &next_set_id, dataTypeNSeq &tot_set) 
#endif
{
	dataTypeNSeq rowIdx;
	dataTypeNSeq colIdx;
	bool rowInSet;
	bool colInSet;
	dataTypeNSeq row_set_id;
	dataTypeNSeq col_set_id;
	
	for (dataTypeNSeq i=0; i<nRows; i++) { //scorri sulle righe calcolate
		rowIdx=firstRow+i;
		row_set_id=reads_to_set[rowIdx]; //0 se non sta in nessun set, k>0 se sta in set k
		if (!row_set_id) { 
			//crea un nuovo insieme con unico elemento l'indice della riga corrente
			tot_set++;
			dataTypeNSeq new_set_id=next_set_id;
			next_set_id++;
			reads_to_set[rowIdx]=new_set_id;
			row_set_id=new_set_id;
		}
		
		//per ogni riga cerca gli elementi maggiori di t
		for (pair<dataTypeNSeq, dataTypeSim> element : simVals[i]) {
			colIdx=element.first; //colIdx può essere >= nReads
				
			#if STEP2
				float normSim; 
				dataTypeSetSim sim=element.second;
				float simf = (float)sim;
				if (sizes[rowIdx] < sizes[colIdx]) normSim=simf/(2*sizes[rowIdx]); //moltiplica per 2 perchè ci sono anche i RC nel DA
				else normSim=simf/(2*sizes[colIdx]);
			#else
				dataTypeSim normSim=element.second;
			#endif
				if (normSim > t) {
					col_set_id=reads_to_set[colIdx%n];
					if (col_set_id && (col_set_id != row_set_id)) { //unisci i due set, SE SONO DIVERSI!!
						if (col_set_id < row_set_id) {
							for (dataTypeNSeq i=0; i<n; i++) {
								if (reads_to_set[i]==row_set_id) reads_to_set[i]=col_set_id;
								else if (reads_to_set[i]>row_set_id) reads_to_set[i]--;
							}
							row_set_id=col_set_id;
							tot_set--;
							next_set_id--;	
						}
						else { //row_set_id < col_set_id
							for (dataTypeNSeq i=0; i<n; i++) {
								if (reads_to_set[i]==col_set_id) reads_to_set[i]=row_set_id;
								else if (reads_to_set[i]>col_set_id) reads_to_set[i]--;
							}
							tot_set--;
							next_set_id--;	
						}
					}
					else if (!col_set_id){ //agiungi al set dove sta l'indice di riga
						reads_to_set[colIdx%n]=row_set_id;
					}
				}
			
		}	
	}
}

int main(int argc, char **argv) {
  
	#if OMP
    double d_total, d_refine;
	#else
    time_t t_refine=0, t_total=0;
    clock_t c_refine=0, c_total=0;
	#endif
	 
	if( argc != 5)
	{
		std::cerr << "Error usage " << argv[0] << " fileFasta t max_rows threads"  << std::endl; // beta non serve più?
		exit(1);
	}

	int num_threads=1;
	sscanf(argv[4], "%d", &num_threads);
	#if OMP
		omp_set_num_threads(num_threads);

    int threads=0;
		#pragma omp parallel
		{
			threads = omp_get_num_threads();
			if(omp_get_thread_num()==0)
				printf("Number of threads: %d\n", threads);
		}
		printf("Number of processors: %d\n", omp_get_num_procs());
	#endif

	string fileFasta=argv[1];
	#if STEP2
		float t;
		sscanf(argv[2], "%f", &t);
	#else
		dataTypeSim t;
		sscanf(argv[2], "%hhu", &t);
	#endif

	dataTypeNSeq max_rows;
	sscanf(argv[3], "%u", &max_rows);
	
	//Read auxiliary file
	dataTypeNSeq numRead;
	dataTypelenSeq minLCP;
	dataTypeNChar maxLen, nClusters;
	bool pos;
	string fileaux=fileFasta.substr(0,fileFasta.find(".fasta"))+".out";
	FILE * outAux = fopen(fileaux.c_str(), "rb");
	if (outAux==NULL) {
		std::cerr << "Error opening " << fileaux << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
	fread(&numRead,sizeof(dataTypeNSeq),1,outAux);
	fread(&minLCP,sizeof(dataTypelenSeq),1,outAux);
	fread(&maxLen,sizeof(dataTypeNChar),1,outAux);
  	fread(&nClusters,sizeof(dataTypeNChar),1,outAux);
  	fread(&pos,sizeof(bool),1,outAux);
	fclose(outAux);
	cout << "numRead: " << numRead << ", minLCP: " << (int)minLCP << ", nClusters: " << nClusters << ", maxLen: " << maxLen << endl;
	//Check maximum cluster length
	if (maxLen>sizeMaxBuf)
	{
    	cerr << "Error Usage: maximum cluster size is " << maxLen << " greater than sizeMaxBuf, please increase sizeMaxBuf in Tools.h" << endl;
		exit(1);
	}
	
	#if STEP2
		//read sets info file, the file contains sets number and size of each set
		string nfInfo=fileFasta.substr(0,fileFasta.find(".fasta"))+".setsinfo";
		FILE * fInfo = fopen(nfInfo.c_str(), "rb");
		if (fInfo==NULL) {
			std::cerr << "Error opening " << nfInfo << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		dataTypeSim old_t;
		fread(&old_t,sizeof(dataTypeSim),1,fInfo);
		dataTypeNSeq numSet;
		fread(&numSet,sizeof(dataTypeNSeq),1,fInfo);
		dataTypeNSeq * sizes = (dataTypeNSeq *)malloc(numSet*sizeof(dataTypeNSeq));
		fread(sizes,sizeof(dataTypeNSeq),numSet,fInfo);
		fclose(fInfo);
		cout << "numSet: " << numSet << endl;
	#endif
	
    //Files .clrs and .da 
	string fnCluster, fnDA;
	#if STEP2
		fnDA=fileFasta+".da.s2";
	#else
		fnDA=fileFasta+".da";
	#endif
	
	std::stringstream ssin;
	ssin << fileFasta << "." << (int)minLCP << ".clrs\0";
	fnCluster=ssin.str();
	
	FILE *InFileCluster[num_threads];
	FILE *InDA[num_threads];
	int t_id=0;
	#if OMP
    for(;t_id<num_threads; t_id++)
    #endif
    {
        InFileCluster[t_id] = fopen(fnCluster.c_str(), "rb");
        if ((InFileCluster[t_id]==NULL)){
			std::cerr << "Error opening " << fnCluster << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		InDA[t_id] = fopen(fnDA.c_str(), "rb");
		if ((InDA[t_id]==NULL)){
			std::cerr << "Error opening " << fnDA << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
	}
	
	#if STEP2
		dataTypeNSeq rem=(numSet)%k;
		dataTypeNSeq q=(numSet)/k;
	#else
		dataTypeNSeq rem=(numRead)%max_rows;
		dataTypeNSeq q=(numRead)/max_rows; 
	#endif
	dataTypeNSeq nIter;
	dataTypeNSeq firstRow=0;
	dataTypeNSeq lastRow=max_rows-1;
	dataTypeNChar AnaCluster;
	
	dataTypeNSeq next_set_id=1;
	dataTypeNSeq tot_set=0;
	#if STEP2
		dataTypeNSeq *reads_to_set = (dataTypeNSeq*) malloc(numSet * sizeof(dataTypeNSeq));
		for (dataTypeNSeq i=0; i<numSet; i++)
			reads_to_set[i]=0;
	#else
		dataTypeNSeq *reads_to_set = (dataTypeNSeq*) malloc(numRead * sizeof(dataTypeNSeq));
		for (dataTypeNSeq i=0; i<numRead; i++)
			reads_to_set[i]=0;
	#endif 
    dataTypeAllSets sets;

	if (rem) 
		nIter=q+1;
	else 
		nIter=q;
	printf("Number of iterations: %d\n", nIter);
	printf("Number of rows in internal memory per iteration: %d\n", max_rows);
	printf("Number of rows in the last iteration (if remainder != 0): %d\n\n", rem);
	
	#if OMP
		d_total = omp_get_wtime();
	#else
		time_start(&t_total, &c_total); //start time
	#endif
	
	#if SQUARE and SAVE_MAT
		string matFileName;
		string matSizesFileName;
		if (pos) {
			matFileName=fileFasta+".mat"+".a"+to_string(minLCP)+"POS.t"+to_string(t);
			matSizesFileName=fileFasta+".matsz"+".a"+to_string(minLCP)+"POS.t"+to_string(t);
		}
		else {
			matFileName=fileFasta+".mat"+".a"+to_string(minLCP)+".t"+to_string(t); 
			matSizesFileName=fileFasta+".matsz"+".a"+to_string(minLCP)+".t"+to_string(t);
		}
		FILE *fpMat;
		FILE *fpMatSizes;
		fpMat=fopen(matFileName.c_str(), "wb");
		if(fpMat==NULL) {
			std::cerr << "Error opening " << matFileName << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		fpMatSizes=fopen(matSizesFileName.c_str(), "wb");
		if(fpMatSizes==NULL) {
			std::cerr << "Error opening " << matSizesFileName << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
	#endif
	
	cerr << "Row analysis started..." << endl;
//------------------------------------------------------------------------------------------------------------- 
	for (dataTypeNSeq i=0; i<nIter; i++) {
		//inizializza il vector di map che contiene i valori di similarità per le righe comprese tra firstRow e lastRow
		dataTypeNSeq nRows=lastRow-firstRow+1;
		#if STEP2
			dataTypeSimsMatS2 simVals;
			for (dataTypeNSeq i=0; i<nRows; i++) {
				simVals.push_back(dataTypeSimsRowS2());
			}
		#else
			dataTypeSimsMat simVals;
			for (dataTypeNSeq i=0; i<nRows; i++) {
				simVals.push_back(dataTypeSimsRow());
			}
		#endif
			
		#if OMP
		  d_refine = omp_get_wtime();
		#else
		  time_start(&t_refine, &c_refine); //start time for clusterRefine
		#endif
	  
		AnaCluster=0;
	  
		//START PARALLEL
		#if OMP
		#pragma omp parallel default(shared) reduction(+:AnaCluster)
		#endif
		{
			int tid = 0;
			int numthreads = 1;
			#if OMP
			tid=omp_get_thread_num();//id_thread
			numthreads = omp_get_num_threads();
			#endif

			dataTypeNChar chunk=(nClusters/numthreads); // ogni thread legge tot clusters
			dataTypeNChar startRead=tid*chunk; // NOTA: startRead qua è diverso da startRead in ClusterLCP.cpp (ogni thread qua iniziava a leggere dall'ultimo elemento del chunk precedente al suo)
			
			dataTypeNChar endRead=(tid+1)*chunk; 
			if(tid==numthreads-1)
				endRead=nClusters;
			
			chunk= endRead-startRead;
			
			fseek(InFileCluster[tid], startRead*sizeof(ElementCluster), SEEK_SET); // ogni thread posiziona il suo puntatore alla parte di file che deve leggere
			
			/**Start analysis for the current rows**/
			#if STEP2
				AnaCluster = clusterAnalyze(chunk,InFileCluster[tid],InDA[tid], numSet, simVals, firstRow, lastRow);
			#else
				AnaCluster = clusterAnalyze(chunk,InFileCluster[tid],InDA[tid], numRead, simVals, firstRow, lastRow);	
			#endif
			/**End analysis for the current rows**/
		}//end-pragma
	  
		#if SMALL 
		// stampa le righe appena calcolate
		cout << "***Row " << firstRow << " to " << lastRow << "***\n" ;
		for (dataTypeNSeq i=0; i<nRows; i++) {
			#if STEP2
			for (pair<dataTypeNSeq, dataTypeSetSim> element : simVals[i]) 
			#else
			for (pair<dataTypeNSeq, dataTypeSim> element : simVals[i]) 
			#endif
			{
				printf("([%d, %d]: %d)", firstRow+i, element.first, (int)element.second);
			}	
			printf("\n");
		}
		#endif
		
		#if SQUARE and SAVE_MAT
			for (dataTypeNSeq i=0; i<nRows; i++) {
				// scrivi size della riga corrente nel file .matsz
				dataTypeNSeq sz=simVals[i].size();
				fwrite(&sz, sizeof(dataTypeNSeq), 1, fpMatSizes);
	
				//scrivi prima il blocco di indici e poi il blocco di valori di similarità relativi alla riga corrente
				dataTypeNSeq * bufferIndeces=(dataTypeNSeq *)malloc(sz*sizeof(dataTypeNSeq));
				dataTypeSim * bufferSims=(dataTypeSim *)malloc(sz*sizeof(dataTypeSim));
				dataTypeNSeq pos=0;
				for (pair<dataTypeNSeq, dataTypeSim> element : simVals[i]) {
					bufferIndeces[pos]=element.first;
					bufferSims[pos]=element.second;
					pos++;
				}
				// i buffer sono pronti, scrivili nel file .mat
				fwrite(bufferIndeces, sizeof(dataTypeNSeq), sz, fpMat); // prima blocco indici
				fwrite(bufferSims, sizeof(dataTypeSim), sz, fpMat); // poi blocco valori di similarità
			}
		#endif
		

		/**Start set analysis for the current rows**/
		#if STEP2
			setAnalysis(simVals, reads_to_set, nRows, firstRow, t, numSet, next_set_id, tot_set, sizes);
		#else
			setAnalysis(simVals, reads_to_set, nRows, firstRow, t, numRead, next_set_id, tot_set);
		#endif
		/**End set analysis for the current rows**/

		#if OMP
			fprintf(stderr,"TIME clusterAnalyze: %.6lf\n", omp_get_wtime()-d_refine);
		#else
			fprintf(stderr,"TIME clusterAnalyze: %.6lf\n", time_stop(t_refine, c_refine));
		#endif

		cout << "Row " << firstRow << " to " << lastRow << " done!" <<"\n" ;
	
		firstRow=lastRow+1;
		if (i==nIter-2 && rem) //sono alla penultima iterazione e c'è il resto
			lastRow=lastRow+rem;
		else  //non sono all'ultima iterazione oppure non c'è il resto
			lastRow=lastRow+max_rows;
	} //end-for
//-------------------------------------------------------------------------------------------------------------  
	#if SQUARE and SAVE_MAT
		fclose(fpMat);
		fclose(fpMatSizes);
	#endif
	
	//metti gli id dei set in [0, totSet-1]
	#if STEP2
	for (dataTypeNSeq i=0; i<numSet; i++) 
	#else
	for (dataTypeNSeq i=0; i<numRead; i++) 
	#endif
	{
		reads_to_set[i]--;
	}
	
	//open output file
	#if STEP2
		string outSets;
		if (pos)
			outSets=fileFasta+".sets2"+".a"+to_string(minLCP)+"POS"+".t"+to_string(old_t)+".t2_"+to_string(t);
		else
			outSets=fileFasta+".sets2"+".a"+to_string(minLCP)+".t"+to_string(old_t)+".t2_"+to_string(t);
	#else
		string outSets;
		if (pos)
			outSets=fileFasta+".sets1"+".a"+to_string(minLCP)+"POS.t"+to_string(t); 
		else
			outSets=fileFasta+".sets1"+".a"+to_string(minLCP)+".t"+to_string(t); 
	
	#endif
	FILE *fpSets;
	string outSetsTxt(outSets.c_str());
	outSetsTxt.append(".txt");
	fpSets=fopen(outSetsTxt.c_str(), "w");
	if(fpSets==NULL) {
		std::cerr << "Error opening " << outSetsTxt << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
	
	//se step 2, leggi il file con i set dello step 1 e scrivi il file con i nuovi set.
	#if STEP2
		//dataTypeNSeq *setsStep1=(dataTypeNSeq *)malloc(numRead*sizeof(dataTypeNSeq));
		FILE * fpsetsStep1;
		string fnsetsStep1;
		if (pos)
			fnsetsStep1=fileFasta+".sets1.a"+to_string(minLCP)+"POS.t"+to_string(old_t)+".txt"; 
		else
			fnsetsStep1=fileFasta+".sets1.a"+to_string(minLCP)+".t"+to_string(old_t)+".txt";
		fpsetsStep1=fopen(fnsetsStep1.c_str(), "r");
		if(fpsetsStep1==NULL) {
			std::cerr << "Error opening " << fnsetsStep1 << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		
		cerr << "Reading " << fnsetsStep1 << endl;
		cerr << "Writing " << outSetsTxt << endl;
		for (dataTypeNSeq i=0; i<numRead; i++) {
			dataTypeNSeq old_set;
			fscanf(fpsetsStep1, "%u", &old_set);
			fprintf(fpSets, "%u %u\n", i, reads_to_set[old_set]);
		}
		fclose(fpSets);
		fclose(fpsetsStep1);
		cerr << "Finished writing " << outSetsTxt << endl;
		cerr << "Finished reading " << fnsetsStep1 << endl;
	#endif 
	
	#if STEP2==0	
		cerr << "Writing " << outSetsTxt << endl;
	#endif
	//l'array con i set è pronto, lo scorri e costruisci la struttura dati per analisi set e, se step 1, scrivi nel file i set_id
	#if STEP2
	for (dataTypeNSeq i=0; i<numSet; i++)
	#else
	for (dataTypeNSeq i=0; i<numRead; i++) 
	#endif
	{	
		#if STEP2==0
			fprintf(fpSets, "%u\n", reads_to_set[i]);
		#endif
		if (sets.count(reads_to_set[i])) { //se il set con questo id è già presente (ci sono già elementi nel set)
			sets[reads_to_set[i]].insert(i);
		}
		else {
			dataTypeSet new_set;
			new_set.insert(i);
			sets.insert({reads_to_set[i], new_set});
		}
	}
	#if STEP2==0
		fclose(fpSets);
		cerr << "Finished writing " << outSetsTxt << endl;
	#endif
	
	//array che contiene in posizione i la size del set i;
	#if STEP2==0
		dataTypeNSeq *sizes = (dataTypeNSeq *)malloc(tot_set*sizeof(dataTypeNSeq));
	#endif
	
	#if SMALL //stampa sets
		printf("***************************************************\n");
		#if STEP2
			printf("t: %f\n", t);
		#else
			printf("t: %d\n", t);
		#endif
		printf("SETS:\n");
		for (pair<dataTypeNSeq, dataTypeSet> s : sets) {
			#if STEP2==0
				sizes[s.first]=s.second.size();
			#endif
			printf("set_id: %d\n", s.first);
			dataTypeSet :: iterator itr;
			for (itr = s.second.begin(); itr != s.second.end(); itr++) 
				cout << (*itr) << "\t";
			cout << endl; 
		}
		printf("Num sets: %d\n", tot_set);
		printf("***************************************************\n");
	#else
		printf("***************************************************\n");
		dataTypeNSeq totR=0;
		dataTypeNSeq size;
		dataTypeNSeq minMaxSize=0;
		dataTypeNSeq nonGrouped=0;
		vector<dataTypeNSeq> maxes;
		vector<dataTypeNSeq>::iterator it;
		
		for (pair<dataTypeNSeq, dataTypeSet> s : sets) {
			#if STEP2
				//scorrere gli elementi di ogni insieme per vedere quali insiemi delle step 1 contiene
				//e sommarne le size (le size le leggi nell'array sizes)
				size=0;
				for (dataTypeSet :: iterator itr = s.second.begin(); itr != s.second.end(); itr++) {
					dataTypeNSeq setSize=sizes[(*itr)];
					size=size+setSize;
				}	
			#else
				size=s.second.size();
				sizes[s.first]=size;
			#endif
			if (size==1) {
				nonGrouped++;
			}
			totR=totR+size;
			if (size>minMaxSize) {
				maxes.push_back(size);
				if (maxes.size()>10) {
					maxes.erase(min_element(maxes.begin(), maxes.end()));
					minMaxSize=*min_element(maxes.begin(), maxes.end());
				}
			}
			else {
				if (maxes.size()<10){
					maxes.push_back(size);
					minMaxSize=size;
				}
			}
		}
		sort(maxes.begin(), maxes.end());
		#if STEP2
			printf("t: %f\n", t);
		#else
			printf("t: %d\n", t);
		#endif
		printf("Num sets: %d, tot reads: %d\n", tot_set, totR);
		printf("First 10 max set size:\n");
		for (int i=maxes.size()-1; i>=0; i--) 
			printf("\t%d\n", maxes[i]);
		printf("Non grouped reads: %d\n", nonGrouped);
		printf("***************************************************\n");
	#endif
	
    //Close .clrs and .da
	t_id=0;
    #if OMP
    for(;t_id<num_threads; t_id++)
    #endif
    {
        fclose(InFileCluster[t_id]);
        fclose(InDA[t_id]);
    }
    #if STEP2==0
		//SCRIVI NUOVO DA-------------------------------------------------------------------------------
		FILE * origDA;
		FILE * newDA;
		
		origDA = fopen(fnDA.c_str(), "rb");
		if (origDA==NULL){
			std::cerr << "Error opening " << fnDA << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		string fnNewDA = fnDA+".s2";
		newDA = fopen(fnNewDA.c_str(), "wb");
		if (newDA==NULL){
			std::cerr << "Error opening " << fnNewDA << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		dataTypeNSeq* DAbuffer= (dataTypeNSeq*)malloc(numRead*sizeof(dataTypeNSeq));
		cerr << "Started writing new DA" << endl;
		dataTypeNChar numcharDA= fread(DAbuffer,sizeof(dataTypeNSeq),numRead,origDA);
		int iter=1;
		while (numcharDA>0) {
			for (dataTypeNChar i=0; i<numcharDA; i++) {
				DAbuffer[i]=reads_to_set[DAbuffer[i]%numRead];
			}
			fwrite(DAbuffer, sizeof(dataTypeNSeq), numcharDA, newDA);
			numcharDA= fread(DAbuffer,sizeof(dataTypeNSeq),numRead,origDA);
			//printf("iter %d done\n", iter);
			//iter++;
		}
		fclose(origDA);
		fclose(newDA);
		free(DAbuffer);
		cerr << "Finished writing new DA" << endl;
		//--------------------------------------------------------------------------------------------- 
		// SCRIVI FILE CON SIZE DEI SET
		//prima scrivi t, poi scrivi numero di set, poi tutte le set size
		cerr << "Started writing sets info file" << endl;
		string fnSetsInfo = fileFasta.substr(0,fileFasta.find(".fasta"))+".setsinfo";
		FILE * setsInfo = fopen(fnSetsInfo.c_str(), "wb");
		if (setsInfo==NULL){
			std::cerr << "Error opening " << fnSetsInfo << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		fwrite(&t,sizeof(dataTypeSim),1,setsInfo);
		//write total number of sets
		fwrite(&tot_set,sizeof(dataTypeNSeq),1,setsInfo);
		//write array with set sizes
		fwrite(sizes, sizeof(dataTypeNSeq),tot_set,setsInfo);
		fclose(setsInfo);
		cerr << "Finished writing sets info file" << endl;
	#endif
    //--------------------------------------------------------------------------------------------------

    cout << "Cluster analysis completed" << "." << endl;
    cout << "Number of clusters: " << AnaCluster << "." << endl;

    #if OMP
      fprintf(stdout,"Time: %.6lf\n", omp_get_wtime()-d_total);
    #else
      fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
    #endif

	return 0;
}

