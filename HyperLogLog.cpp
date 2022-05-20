#include "HyperLogLog.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include "wyhash32.h"

using namespace std;
#define seed 5 //seed para hash
#define lim 4294967296 //2^32

HyperLogLog::HyperLogLog(unsigned char n1, unsigned char n2,unsigned char k,int numThreads, int n){
	omp_set_num_threads(numThreads);
	p=n1;
	b=n2;
	bits_v2=(1<<b)-1; //2^b -1
	N=1<<p; //2^V1
	a_m=(0.7213/(1+(1.079/N)))*N*N;
	for(unsigned char i=0;i<12;++i)
		bit_mask.emplace_back(~((ullint)0x1F<<(5*i)));
	cards.resize(n);
	sketch.resize(n);
	names.resize(n);
	kmer_length=to_string((int)k);
	sketch_size=to_string(p);
}
HyperLogLog::~HyperLogLog(){}

void HyperLogLog::addSketch(string genome,int i){
	//inicializa sketch en 0s, y añade su nombre a un vector
	for(unsigned int j=0;j<N/12;++j) //cada celda tendra 12 buckets del sketch array
		sketch[i].emplace_back(0);
	if((long)N%12)
		sketch[i].emplace_back(0);
	names[i]=genome;
}

void HyperLogLog::insert(ullint kmer,int i){
	//calculamos el hash usando wyhash de 32 bits (hashing muy rapido)
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	unsigned int v1=hash>>b; //se desplaza el hash b bits a la derecha
	ullint v2=hash&bits_v2; // AND 2^b -1
	if(v2==0) v2=b+1;
	else v2=__builtin_clz(v2)+1-p; //cuenta cantidad de bits 0 mas significativo
	//lo cual sirve para encontrar posicion de 1 mas a la izquierda
	//se resta p, ya que al trabajar con int, se tiene 32 bits,
	//y siempre habran al menos p ceros a la izquierda
	unsigned int indice=v1/12; //indice de celda
	unsigned char bucket=v1%12; //12 buckets por celda
	ullint temp=(sketch[i][indice]>>(5*bucket))&0x1F;
	if(v2 > temp) sketch[i][indice]=(sketch[i][indice]&(bit_mask[bucket]))|(v2<<(5*bucket));
}

void HyperLogLog::loadSketch(char *filename){ //¿que pasa si tamaño del sketch es distinto al de argumento?
	/*FILE *fp = fopen(filename,"r");
	string genome=filename; //se debe determinar nombre del genoma a partir del archivo .hll
	genome.erase(genome.find(".w."),string::npos);
	vector<ullint> sketch;
	unsigned int ceros=N;
	while(!feof(fp)){
		ullint celda;
		int ret=fscanf(fp,"%llu",&celda);
		sketch.emplace_back(celda);
	}
	fclose(fp);
	sketches.push_back(pair<string,vector<ullint>>(genome,sketch));*/
}

void HyperLogLog::saveSketches(){
	int tam=sketch.size();
	#pragma omp parallel for
	for(int i=0;i<tam;++i){
		vector<ullint> s=sketch[i];
		string temp=names[i];
		temp+=".w."+kmer_length+".spacing."+sketch_size+".hll";
		char* filename=(char*)temp.c_str();
		FILE *fp = fopen(filename,"w");
		for(vector<ullint>::iterator it=s.begin();it!=s.end();it++)
			fprintf(fp,"%llu\n",*it);
		fclose(fp);
	}
}

void HyperLogLog::estCard(){
	int tam=sketch.size();
	#pragma omp parallel for
	for (int j = 0; j < tam; ++j){
		int ceros=N;
		vector<unsigned int> w;
		for(unsigned char i=0;i<b+2;++i)
			w.emplace_back(0);
		vector<ullint>::iterator it1=sketch[j].begin();
		vector<ullint>::iterator fin=sketch[j].end();
		//contamos las repeticiones de w en cada sketch
		while(it1!=fin-1){ 
			ullint i1=*it1;
			for(char i=0;i<12;++i){ //12 registros por celda
				ullint temp1=i1&0x1F;
				if(temp1) ceros--;
				w[temp1]++;
				i1=i1>>5;
			}
			++it1;
		}
		unsigned char tam=(long)N%12;
		ullint i1=*it1;
		for(char i=0;i<tam;++i){ //12 registros por celda
			ullint temp1=i1&0x1F;
			w[temp1]++;
			i1=i1>>5;
		}
		//eliminamos las repeticiones por celda extra vacía creada (?)
		w[0]-=12;
		long double card=0;
		for(unsigned char i=0;i<b+2;++i)
			if(w[i]) card+=(long double)w[i]/(long double)(1<<i);
		//media armonica
		card=(long double)a_m/card;
		if(ceros && card<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
			//printf("ceros A:%u\n",cerosA);
			printf("linear counting\n");
			card=N*log(N/ceros);
		}
		else if(card>lim/30){
			printf("valor muy grande\n");
			card=-lim*log(1-(card/lim));
		}
		printf("estimacion cardinalidad %s: %Lf\n",names[j].c_str(),card);
		cards[j]=card;
	}
}

void HyperLogLog::estJaccard(){ //recibe numero de genomas
	int tam=sketch.size();
	#pragma omp parallel for //no se puede usar collapse aqui
	for(int j1=0;j1<tam;j1++){
		for(int j2=j1+1;j2<tam;j2++){
			vector<ullint>::iterator it1=sketch[j1].begin();
			vector<ullint>::iterator it2=sketch[j2].begin();
			vector<ullint>::iterator fin=sketch[j1].end();
			int ceros=N;
			vector<unsigned int> wU;
			for(unsigned char i=0;i<b+2;++i)
				wU.emplace_back(0);
			//contamos las repeticiones de w en cada sketch
			while(it1!=fin-1){ 
				ullint i1=*it1,i2=*it2;
				for(char i=0;i<12;++i){ //12 registros por celda
					ullint temp1=i1&0x1F,temp2=i2&0x1F;
					if(temp1||temp2) --ceros;
					(temp1>temp2) ? wU[temp1]++ : wU[temp2]++;
					i1=i1>>5;
					i2=i2>>5;
				}
				++it1;
				++it2; //avanza en tabla B
			}
			unsigned char tam=(long)N%12;
			ullint i1=*it1,i2=*it2;
			for(char i=0;i<tam;++i){ //12 registros por celda
				ullint temp1=i1&0x1F,temp2=i2&0x1F;
				if(temp1||temp2) --ceros;
				(temp1>temp2) ? wU[temp1]++ : wU[temp2]++;
				i1=i1>>5;
				i2=i2>>5;
			}
			//eliminamos las repeticiones por celda extra vacía creada (?)
			wU[0]-=12;
			
			long double cardU=0;
			for(unsigned char i=0;i<b+2;++i){
				if(wU[i]) cardU+=(long double)wU[i]/(long double)(1<<i);
			}
			//media armonica
			cardU=(long double)a_m/cardU;
			if(ceros && cardU<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
				//printf("ceros A U B:%u\n",cerosU);
				printf("linear counting\n");
				cardU=N*log(N/ceros);
			}
			else if(cardU>lim/30){
				printf("valor muy grande\n");
				cardU=-lim*log(1-(cardU/lim));
			}
			//printf("estimacion cardinalidad %s U %s: %Lf\n",names[j1].c_str(),names[j2].c_str(),cardU);

			long double jaccard=(cards[j1]+cards[j2]-cardU)/cardU;
			if(jaccard<0) jaccard=0;
			printf("estimacion jaccard %s U %s: %Lf\n",names[j1].c_str(),names[j2].c_str(),jaccard);
		}
	}
}
