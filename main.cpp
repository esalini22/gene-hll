#include <stdio.h>
#include <string.h>
#include <string>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <omp.h>
#include <thread>
#include <vector>
#include <algorithm>
#include <vector>
#include <immintrin.h> 
#include "HyperLogLog.h"
#define lim 4294967296 //2^32
#define lim_int 2147483647

using namespace std;
typedef unsigned long long int ullint;

unsigned char k=31; //largo de kmer
ullint bits_G;
ullint bits_T;
ullint bits_C;
ullint BITS;

inline float hsum_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 maxs = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, maxs); // high half -> low half
    maxs        = _mm_add_ss(maxs, shuf);
    return        _mm_cvtss_f32(maxs);
}

inline float hsum_avx(__m256 v) {
    __m128 lo = _mm256_castps256_ps128(v);   // low 128
    __m128 hi = _mm256_extractf128_ps(v, 1); // high 128
           lo = _mm_add_ps(lo, hi);          // max the low 128
    return hsum_sse3(lo);                    // and inline the sse3 version
}

vector<vector<float>> estJaccard(vector<HyperLogLog*> hll, vector<float> &cards,int b,int N, int numThreads){
	omp_set_num_threads(numThreads);
	__m256i vec3; //vector de mascaras de bits
	ullint bits_and[4]={15,15,15,15};
	//vec3=_mm256_loadu_si256((const __m256i*)&bits_and[0]);

	float a_m=(0.7213/(1+(1.079/N)))*N*N;
	int ciclos_red=(b+2)/8+(((b+2)%8)>0);

	int tam=hll.size();
	vector<vector<float>> jaccards(tam);

	#pragma omp parallel for //no se puede usar collapse aqui
	for(int j1=0;j1<tam;j1++){
		vector<float> temp_jaccards; //para evitar false sharing
		temp_jaccards.resize(tam-j1-1);
		for(int j2=j1+1;j2<tam;j2++){
			vector<ullint> s1ref=(hll[j1]->getSketch());
			vector<ullint> s2ref=(hll[j2]->getSketch());
			vector<ullint>::iterator it1=s1ref.begin();
			vector<ullint>::iterator it2=s2ref.begin();
			vector<ullint>::iterator fin=s1ref.end();
			vector<float> wU(32,0.0);
			//contamos las repeticiones de w en cada sketch
			while(it1!=fin){ 
				ullint i1=*it1,i2=*it2;

				//solo funciona con avx 512, por lineas _mm256_max_epu64()
				/*ullint it_array1[16]={i1,(i1>>4),(i1>>8),(i1>>12),(i1>>16),(i1>>20),(i1>>24),(i1>>28),(i1>>32),(i1>>36),(i1>>40),(i1>>44),(i1>>48),(i1>>52),(i1>>56),(i1>>60)};
				ullint it_array2[16]={i2,(i2>>4),(i2>>8),(i2>>12),(i2>>16),(i2>>20),(i2>>24),(i2>>28),(i2>>32),(i2>>36),(i2>>40),(i2>>44),(i2>>48),(i2>>52),(i2>>56),(i2>>60)};
				__m256i vec4[4],vec5[4];
				vec4[0]=_mm256_loadu_si256((const __m256i*)&it_array1[0]);
				vec4[1]=_mm256_loadu_si256((const __m256i*)&it_array1[4]);
				vec4[2]=_mm256_loadu_si256((const __m256i*)&it_array1[8]);
				vec4[3]=_mm256_loadu_si256((const __m256i*)&it_array1[12]);
				vec4[0]=_mm256_and_si256(vec3,vec4[0]);
				vec4[1]=_mm256_and_si256(vec3,vec4[1]);
				vec4[2]=_mm256_and_si256(vec3,vec4[2]);
				vec4[3]=_mm256_and_si256(vec3,vec4[3]);
				vec5[0]=_mm256_loadu_si256((const __m256i*)&it_array2[0]);
				vec5[1]=_mm256_loadu_si256((const __m256i*)&it_array2[4]);
				vec5[2]=_mm256_loadu_si256((const __m256i*)&it_array2[8]);
				vec5[3]=_mm256_loadu_si256((const __m256i*)&it_array2[12]);
				vec5[0]=_mm256_and_si256(vec3,vec5[0]);
				vec5[1]=_mm256_and_si256(vec3,vec5[1]);
				vec5[2]=_mm256_and_si256(vec3,vec5[2]);
				vec5[3]=_mm256_and_si256(vec3,vec5[3]);

				vec4[0]=_mm256_max_epu64(vec4[0],vec5[0]);
				vec4[1]=_mm256_max_epu64(vec4[1],vec5[1]);
				vec4[2]=_mm256_max_epu64(vec4[2],vec5[2]);
				vec4[3]=_mm256_max_epu64(vec4[2],vec5[3]);
				__attribute__ ((aligned (32))) ullint out[4];
				for(char c=0;c<4;++c){
					_mm256_store_si256((__m256i *)&out[0],vec4[c]);
					if(out[0]) ceros--;
					wU[out[0]]++;
					if(out[1]) ceros--;
					wU[out[1]]++;
					if(out[2]) ceros--;
					wU[out[2]]++;
					if(out[3]) ceros--;
					wU[out[3]]++;
				}*/

				for(char i=0;i<16;++i){ //16 registros por celda
					ullint temp1=i1&0xF,temp2=i2&0xF;
					(temp1>temp2) ? wU[temp1]++ : wU[temp2]++;
					i1=i1>>4;
					i2=i2>>4;
				}
				++it1;
				++it2; //avanza en tabla B
			}
			
			float cardU=0.0;

			//for(unsigned char i=0;i<b+2;++i)
			//	if(wU[i]) cardU+=(float)wU[i]/(float)(1<<i);

			float w2[32];
			for(int i=0;i<32;++i) w2[i]=1.0;
			int respow=1;
			for(int i=0;i<b+2;++i){
				w2[i]=(float)respow;
				respow=respow<<1;
			}
			for(int i=0;i<32;++i){
				printf("%d. %f %f\n",i,wU[i],w2[i]);
			}

			__m256 vec,vec2;
			for(int i=0;i<ciclos_red;++i){
				vec=_mm256_loadu_ps((const float *)&wU[i*8]);
				vec2=_mm256_loadu_ps((const float *)&w2[i*8]);
				vec=_mm256_div_ps(vec,vec2);
				cardU+=hsum_avx(vec);
				printf("sum: %f\n",hsum_avx(vec));
			}

			int ceros = wU[0];

			//media armonica
			cardU=(float)a_m/cardU;
			if(ceros && cardU<=5*N/2) //C_HLL, ln cuando hay muchos ceros;
				cardU=N*log(N/ceros);
			else if(cardU>lim/30)
				cardU=-lim*log(1-(cardU/lim));
			printf("estimacion cardinalidad union: %f ceros: %d\n",cardU,ceros);

			float jaccard=(cards[j1]+cards[j2]-cardU)/cardU;
			if(jaccard<0) jaccard=0;
			temp_jaccards[j2-j1-1]=jaccard;
		}
		jaccards[j1]=temp_jaccards;
	}
	return jaccards;
}

void printMatrix(vector<vector<float>> &jaccards, vector<string> names){
	int tam=names.size();
	//printf("	");
	char guion='-';
	for(int j=0;j<tam;j++)
		printf("%3s ",names[j].c_str());
	printf("\n");
	for(int j1=0;j1<tam;j1++){
		printf("%3s ",names[j1].c_str());
		for(int j2=0;j2<tam;j2++){
			if(j2<j1+1){
				//printf("- ");
				printf("%3c ",guion);
				continue;
			}
			printf("%3f ",jaccards[j1][j2-j1-1]);
		}
		printf("\n");
	}
}

void saveOutput(char* filename,vector<string> names, vector<vector<float>> jaccards){ //guarda la matriz en txt
	int tam=names.size();
	FILE *fp=fopen(filename,"w");
	fprintf(fp,"	");
	for(int j=0;j<tam;j++)
		fprintf(fp,"%s ",names[j].c_str());
	fprintf(fp,"\n");
	for(int j1=0;j1<tam;j1++){
		fprintf(fp,"%s ",names[j1].c_str());
		for(int j2=0;j2<tam;j2++){
			if(j2<j1+1){
				fprintf(fp,"- ");
				continue;
			}
			fprintf(fp,"%f ",jaccards[j1][j2-j1-1]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


/*void to_kmer(unsigned long long int num,unsigned char k){
	unsigned long long int temp=num;
	char kmer[k+1];
	kmer[k]='\0';
	for(int i=k-1;i>=0;--i){
		char c=temp&3;
		if(c==0) kmer[i]='A';
		else if(c==1) kmer[i]='C';
		else if(c==2) kmer[i]='G';
		else if(c==3) kmer[i]='T';
		temp=temp>>2;
	}
	printf("num: %llx\n",num);
	printf("kmer: %s\n",kmer);
}*/

void leer(char *genome,HyperLogLog *hll){
	char c; //para lectura
	ullint kmer=0,comp=0;
	string linea;
	ifstream indata(genome);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",genome);
		exit(1);
	}

	hll->addSketch(genome); //añade el sketch inicializado con 0s

	getline(indata,linea); //salta primera linea
	getline(indata,linea); //para primer kmer
	string::iterator it=linea.begin();

	//al leerse cada caracter, se insert la base al final del kmer no canonico
	//y luego se insert su complemento al inicio del complemento del reverso del kmer
	//el inicio estara dado por el largo del kmer
	for(unsigned char j=0;j<k;++j){ //inicializacion - primer kmer
		kmer=kmer<<2;
		comp=comp>>2;
		c=*it;
		++it;
		if(c=='A') comp=comp|bits_T; //A=00 0xC000000000
		else if(c=='C'){ //C=01
			kmer=kmer|0x1;
			comp=comp|bits_G; //0x8000000000
		}
		else if(c=='G'){ //G=10
			kmer=kmer|0x2;
			comp=comp|bits_C; //0x4000000000
		}
		else if(c=='T') kmer=kmer|0x3; //T=11
	}
	(kmer>comp) ? hll->insert(comp) : hll->insert(kmer); //insert kmer canonico

	while(!indata.eof()){
		while(it!=linea.end()){
			c=*it;
			if(c=='A' || c=='C' || c=='G' || c=='T'){ //lee nueva base
				kmer=(kmer<<2)&BITS; //desplaza bases a la izquierda
				comp=comp>>2; //desplaza bases a la derecha
				if(c=='A') comp=comp|bits_T; //A=00
				else if(c=='C'){ //C=01
					kmer=kmer|0x1;
					comp=comp|bits_G;
				}
				else if(c=='G'){ //G=10
					kmer=kmer|0x2;
					comp=comp|bits_C;
				}
				else if(c=='T') kmer=kmer|0x3; //T=11

				(kmer>comp) ? hll->insert(comp) : hll->insert(kmer); //insert kmer canonico
			}
			else if(c=='>') break;
			++it;
		}
		getline(indata,linea);
		it=linea.begin();
	}
	indata.close();
}

//obtiene archivos de un txt
//formato de archivos:
//genoma1
//genoma2
//etc. (1 genoma por linea)
vector<string> readFromFile(char* paths){
	ifstream indata(paths);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",paths);
		exit(1);
	}
	vector<string> genomes;
	string filename;
	while(!indata.eof()){
		getline(indata,filename);
		if(filename!="") genomes.push_back(filename);
	}
	indata.close();
	return genomes;
}

//obtiene archivos de la linea de argumentos
vector<string> getPaths(char** argv, int argc){
	vector<string> genomes;
	for(int i=1;i<argc;++i){
		if(!strcmp(argv[i],"-k") || !strcmp(argv[i],"-p") || !strcmp(argv[i],"-t") || !strcmp(argv[i],"-o") || !strcmp(argv[i],"-d") || !strcmp(argv[i],"-r")) ++i;
		else if(strcmp(argv[i],"-s")) genomes.push_back(argv[i]);
	}
	return genomes;
}

vector<string> readCompressedFromFile(char* paths){
	ifstream indata(paths);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",paths);
		exit(1);
	}
	vector<string> genomes;
	string filename;
	while(!indata.eof()){
		getline(indata,filename);
		if(filename!="") genomes.push_back(filename);
	}
	indata.close();
	return genomes;
}

//obtiene archivos de la linea de argumentos
vector<string> getCompressed(char** argv, int argc){
	vector<string> genomes;
	for(int i=1;i<argc;++i){
		if(!strcmp(argv[i],"-k") || !strcmp(argv[i],"-p") || !strcmp(argv[i],"-t") || !strcmp(argv[i],"-o") || !strcmp(argv[i],"-f") || !strcmp(argv[i],"-r")) ++i;
		else if(!strcmp(argv[i],"-d")) genomes.push_back(argv[i+1]);
	}
	return genomes;
}

//formato: ./hll -opcion valor genomas
//o bien ./hll genomas -opcion valor
//no detecta caso en que se introduza opcion o valor invalido
int main(int argc, char *argv[]){
	if(argc<3) {
		printf("No hay suficientes argumentos\n");
		exit(1);
	}
	unsigned char p=12;
	char** option;
	char** end=argv+argc;
	option=std::find((char**)argv,end,(const std::string&)"-k"); //cambia largo de kmer
	if(option!=end){
		char val=atoi(*(option+1));
		if(val<32 && val>19) k=val;
	}
	option=std::find((char**)argv,end,(const std::string&)"-p"); //cambia tamaño de sketch
	if(option!=end){
		char val=atoi(*(option+1));
		if(val<32 && val>8) p=val;
	}
	vector<string> genomes,compressed;
	option=std::find((char**)argv,end,(const std::string&)"-f"); //lee de txt
	if(option!=end) genomes=readFromFile(*(option+1));
	else genomes=getPaths(argv,argc);

	option=std::find((char**)argv,end,(const std::string&)"-r"); //lee de txt
	if(option!=end) compressed=readCompressedFromFile(*(option+1));
	else compressed=getCompressed(argv,argc);

	int tam=genomes.size(),tam2=compressed.size();
	printf("tam: %d, tam2: %d\n",tam,tam2);

	for(vector<string>::iterator it=genomes.begin();it!=genomes.end();++it)
		printf("%s ",(char*)(*it).c_str());

	printf("\nk: %d p: %d\n",k,p);
	

	int numThreads=min(tam+tam2,(int)std::thread::hardware_concurrency());
	printf("numThreads: %d, tam: %d, maxThreads: %d\n",numThreads,tam+tam2,(int)std::thread::hardware_concurrency());

	option=std::find((char**)argv,end,(const std::string&)"-t"); //asigna manualmente la cantidad de hebras
	if(option!=end) numThreads=atoi(*(option+1));

	printf("threads: %d\n",numThreads);

	vector<HyperLogLog*> v_hll;
	for(int i=0;i<tam+tam2;++i){
		HyperLogLog *hll;
		hll = new HyperLogLog(p,32-p,k);
		v_hll.push_back(hll);
	}

	//se trabaja a nivel de bits, es mas rapido que trabajar con strings
	//aca se determina como se insertaran las bases complementarias en el complemento del reverso del kmer
	//es decir, se determina donde esta el inicio (en bits) de dicho kmer
	const ullint desp=(2*(k-1));
	bits_G=(ullint)2<<desp;
	bits_T=(ullint)3<<desp;
	bits_C=(ullint)1<<desp;
	//esto sirve para eliminar la primera base del kmer, luego de desplazar las bases a la izquierda para leer la nueva base
	BITS=(bits_C-1)<<2;

	//lee paralelamente cada archivo
	omp_set_num_threads(numThreads);
	#pragma omp parallel
	{
		#pragma omp single
		{
			for(int i=0;i<tam;i++){
				#pragma omp task
				leer((char*)genomes[i].c_str(),v_hll[i]);
			}
			for(int i=0;i<tam2;i++){
				#pragma omp task
				v_hll[i+tam]->loadSketch((char*)compressed[i].c_str());
			}
		}
	}
	
	option=std::find((char**)argv,end,(const std::string&)"-s"); //guarda los sketches
	if(option!=end){
		for(int i=0;i<tam+tam2;++i)
			v_hll[i]->saveSketch();
	}
	vector<string> names(tam+tam2);
	vector<float> cards(tam+tam2);
	#pragma omp parallel for
	for(int i=0;i<tam+tam2;i++){
		//printf("%d\n",i);
		names[i]=v_hll[i]->getName();
		cards[i]=v_hll[i]->estCard();
		printf("getcard\n");
	}
	printf("jaccard\n");
	vector<vector<float>> jaccards=estJaccard(v_hll,cards,32-p,1<<p,numThreads);
	/*option=std::find((char**)argv,end,(const std::string&)"-o"); //guarda la matriz en txt
	if(option!=end) saveOutput(*(option+1));
	else */printMatrix(jaccards,names);

	for(int i=0;i<tam+tam2;++i)
		delete v_hll[i];

	return 0;
}
