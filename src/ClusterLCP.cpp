#include "Tools.h"

#ifndef OMP
	#define OMP 1
#endif

#ifndef POS
	#define POS 0
#endif

/*
Step 1. detects alpha-clusters, i.e., blocks containing symbols belonging both to reads and to genomes, and whose associated suffixes share a common context of minimum length alpha.

As preprocessing, it need to have the two data structures fileFasta.lcp and fileFasta.da computed.

Input: fileFasta, total number of reads, total number of genomes, alpha, threads.
		
Output: fileFasta.alpha.clrs containing a pair ElementCluster (pStart,len) for each alpha-cluster detected; one auxiliary file.
*/

#if POS
void StartCluster(dataTypeNChar ind, ElementCluster &cluster)
{
	cluster.pStart = ind-1;
}
#endif

void StartOrRemain(dataTypeNChar ind,bool &init, ElementCluster &cluster)
{	
	if (not init) //Start a new cluster
	{
		init=true;
		cluster.pStart = ind-1;
	}
}

int Close(ElementCluster &cluster, dataTypeNChar ind, dataTypeNChar &lengthClust,vector<ElementCluster> &vOutput)
{
	cluster.len = ind - cluster.pStart;
	if (lengthClust<cluster.len)
		lengthClust=cluster.len;
        
	vOutput.push_back(cluster);
  
	return 1;
}
                          


int main(int argc, char **argv) {

	#if OMP
		double d_total;
	#else
		time_t t_refine=0, t_total=0;
		clock_t c_refine=0, c_total=0;
	#endif

	if( argc != 5 )
	{	
		std::cerr << "Error usage: " << argv[0] << " fileFasta numReads alpha threads" << std::endl; 
		exit(1);
	}
    
	string fileFasta=argv[1]; 
	dataTypelenSeq alpha;
	dataTypeNSeq numReads;
	
	sscanf(argv[2], "%u", &numReads);
	sscanf(argv[3], "%hhu", &alpha); // se dataTypelenSeq e' uchar %hhu, se uint %u
	
    int num_threads=1;
    #if OMP
		sscanf(argv[4], "%d", &num_threads);
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
	
	//Open files
	string fnLCP, fileOutput; 
	fnLCP = fileFasta+".lcp\0";
	std::stringstream ssout;
	ssout << fileFasta << "." << (int)alpha << ".clrs\0";
	fileOutput=ssout.str();
	
	FILE *OutCluster=fopen(fileOutput.c_str(), "w");
	if(OutCluster==NULL) {
		cerr << "Error opening " << fileOutput << ".";
		exit(1);
	}
	
	FILE *InLCP[num_threads];
	int tid=0;
	#if OMP
    for(;tid<num_threads; tid++)
    #endif
    {
        InLCP[tid] = fopen(fnLCP.c_str(), "rb");
        if(!tid){
			std::cout << "\n\t" << fnLCP;
			if ((InLCP[tid]==NULL)){
				std::cerr << "Error opening " << fnLCP << "." << std::endl;
				printf("fopen failed, errno = %d\n", errno);
				exit (EXIT_FAILURE);
			}
        }
    }
	
	//InLCP dimension
	fseek(InLCP[0], 0, SEEK_END);
	dataTypeNChar sizeInLCP=ftell(InLCP[0])/sizeof(dataTypelenSeq);
	//Clustering
	#if OMP
		d_total = omp_get_wtime();
	#else
		time_start(&t_total, &c_total); //start time
	#endif
	
	dataTypeNChar maxLen=0, nClusters=0;
	
	//START PARALLEL
    #if OMP
    #pragma omp parallel default(shared) reduction(max:maxLen) reduction(+:nClusters)
    #endif
    {
        int t_id = 0;
		int numthreads = 1;
    #if OMP
        t_id=omp_get_thread_num();//id_thread
        numthreads = omp_get_num_threads();
    #endif
        
        dataTypeNChar chunk=(sizeInLCP/numthreads);
		dataTypeNChar bufferSize=(BUFFERLCPSIZE/numthreads);
		assert(bufferSize>0);
        
        dataTypeNChar startRead=t_id*chunk-1;
        dataTypeNChar endRead=(t_id+1)*chunk;
		
		if(t_id==0)
			startRead=0;
		
        if(t_id==numthreads-1)
            endRead=sizeInLCP;
    #if OMP
		double start=omp_get_wtime();
    #endif
        
		fseek(InLCP[t_id], startRead*sizeof(dataTypelenSeq), SEEK_SET);

		//Output file contains a collection of pairs (pStart, len) each one corresponding to one alpha-cluster
		ElementCluster cluster;
		cluster.pStart=0;
		cluster.len=0;
		vector<ElementCluster> vOutput;
        
		bool init=false; //init is true if a cluster is open
		
	#if POS
		dataTypelenSeq pred = -1;
		bool grow = true;
		dataTypeNChar cont = 0;
	#endif
	
		//To read LCP file
		dataTypeNChar numcharLCP;
		dataTypelenSeq* bufferLCP = new dataTypelenSeq[bufferSize];
		
		dataTypeNChar index=1; 
		if(t_id!=0)
			index=t_id*chunk; //startRead-1
		
		numcharLCP=fread(&bufferLCP[0],sizeof(dataTypelenSeq),1,InLCP[t_id]);

		chunk=endRead-startRead-1; //per il primo thread vale 1 in meno rispetto agli altri threads. Indica quanti valori devo ancora leggere
		
		while ((chunk>0) && (bufferLCP[0]>=alpha))
		{
			numcharLCP=fread(&bufferLCP[0],sizeof(dataTypelenSeq),1,InLCP[t_id]);
			chunk--;
			index++;
		}
		
		if(bufferLCP[0]<alpha) 
		{
			while(numcharLCP>0)
			{
				if(bufferSize>=chunk)
					bufferSize=chunk;
					
				numcharLCP = fread(bufferLCP,sizeof(dataTypelenSeq),bufferSize,InLCP[t_id]);
                
				for(dataTypeNChar indexbuffer=0; indexbuffer<numcharLCP; indexbuffer++)
				{
				#if POS	
					if (bufferLCP[indexbuffer]>=alpha) {
						if (not init) {
							StartCluster(index, cluster); 
							init=true;
							grow=true;
							cont=0;
						}
						else { //cluster già iniziato, o ci resti dentro o lo chiudi perchè ti accorgi di aver trovato un minimo e apri subito un altro cluster
							if (grow && bufferLCP[indexbuffer]<pred) grow = false;
							else if ((not grow) && bufferLCP[indexbuffer]<pred) cont=0;
							else if ((not grow) && bufferLCP[indexbuffer]==pred) cont++;
							else if ((not grow) && bufferLCP[indexbuffer]>pred) {
								nClusters+=Close(cluster, index-1-cont, maxLen,vOutput);
								init=false;
								cluster.pStart=0, cluster.len=0;
								StartCluster(index-cont, cluster);
								init=true;
								grow=true; 
								cont=0;
							}
						}
						pred=bufferLCP[indexbuffer];
					}
				#else
					if(bufferLCP[indexbuffer]>=alpha) //Start or remain in a cluster
						StartOrRemain(index, init, cluster);
				#endif	
					else    //End a cluster
					{
						if (init)
							nClusters+=Close(cluster, index, maxLen,vOutput);
						
						init=false;
						cluster.pStart=0, cluster.len=0;
					}
					index++;
				}//end-for
				
                #if OMP
					#pragma omp critical
				#endif
                {
                    for(dataTypeNChar i=0;i<vOutput.size();i++) {
                        fwrite(&vOutput[i],sizeof(ElementCluster),1,OutCluster);
						//printf("start: %d, len: %d\n", (int)vOutput[i].pStart, (int)vOutput[i].len);
                    }
                }
                
				//Read the LCP file
				chunk-=numcharLCP;
                vOutput.clear();
			}//end-while
			
			//We need to close a possibly open cluster 
			if(init && (t_id==num_threads-1)) //Last thread 
				nClusters+=Close(cluster, index, maxLen,vOutput);
			else if (init && (t_id!=num_threads-1)) //Manage straddling clusters 
			{	
				bool stayIn=true;
				do {
					numcharLCP=fread(&bufferLCP[0],sizeof(dataTypelenSeq),1,InLCP[t_id]);
					if(numcharLCP>0 && bufferLCP[0]>=alpha) {//stai nel cluster o chiudi perchè hai trovato un minimo
						//se no POS non fai nulla, se POS simile a sopra ...
					#if POS
						if (grow && bufferLCP[0]<pred) grow = false;
						else if ((not grow) && bufferLCP[0]<pred) cont=0;
						else if ((not grow) && bufferLCP[0]==pred) cont++;
						else if ((not grow) && bufferLCP[0]>pred) {
							nClusters+=Close(cluster, index-1-cont, maxLen,vOutput);
							init=false;
							cluster.pStart=0, cluster.len=0;
							StartCluster(index-cont, cluster);
							init=true;
							grow=true;
							cont=0;
						}
						pred=bufferLCP[0];
					#endif
					
					}
					else//End THE cluster
					{
						stayIn=false;
						nClusters+=Close(cluster, index, maxLen,vOutput);
					}//end-else
					index++;
				} while(stayIn);
				
			}//end-else-if
            #if OMP
				#pragma omp critical
			#endif
			
            {
                for(dataTypeNChar i=0;i<vOutput.size();i++){
                    fwrite(&vOutput[i],sizeof(ElementCluster),1,OutCluster);
					//printf("start: %d, len: %d\n", (int)vOutput[i].pStart, (int)vOutput[i].len);
                }
            }
		}//end-if
		#if OMP
			#pragma omp critical
			{
				std::cerr << "TIME THREAD " << t_id << " = " << omp_get_wtime()-start << "(in seconds)\n";
			}
		#endif
		
		delete[] bufferLCP;
		
	}//end-pragma
	tid=0;
    #if OMP
    for(;tid<num_threads; tid++)
    #endif
    {
        fclose(InLCP[tid]);
    }
	fclose(OutCluster);
	
	string fileaux=fileFasta.substr(0,fileFasta.find(".fasta"))+".out";
	
   //Write auxiliary file
	FILE * outAux = fopen(fileaux.c_str(), "w");
	if (outAux==NULL) {
		std::cerr << "Error opening " << fileaux << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
  
	fwrite(&numReads,sizeof(dataTypeNSeq),1,outAux);
	fwrite(&alpha,sizeof(dataTypelenSeq),1,outAux);
	fwrite(&maxLen,sizeof(dataTypeNChar),1,outAux);
	fwrite(&nClusters,sizeof(dataTypeNChar),1,outAux);
	#if POS
		bool pos = true;
		fwrite(&pos,sizeof(bool),1,outAux);
	#else
		bool pos = false;
		fwrite(&pos,sizeof(bool),1,outAux);
	#endif
                                  
	fclose(outAux);
                                  
	cout << "Clustering process with alpha=" << (int)alpha << " completed.\nTotal number of clusters: " << nClusters << ".\nMaximum cluster size: " << maxLen << "." << endl;
	#if POS
		cout << "Positional Clustering" << endl;
	#endif
	#if OMP
		fprintf(stdout,"Time: %.6lf\n", omp_get_wtime()-d_total);
    #else
		fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
    #endif

return 0;
}

