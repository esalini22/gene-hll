include <stdio.h>
#include <string.h>
#include <string>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <omp.h>
#include <thread>
#include <vector>
#include <algorithm>
#include "HyperLogLog.h"
using namespace std;
typedef unsigned long long int ullint;

HyperLogLog *hll;
unsigned char k=31; //largo de kmer
ullint bits_G;
ullint bits_T;
ullint bits_C;
ullint BITS;

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

void leer(char *genome,int i){
	char c; //para lectura
	ullint kmer=0,comp=0;
	string linea;
	ifstream indata(genome);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",genome);
		exit(1);
	}

	/*vector<ullint>::iterator sketch=*/hll->addSketch(genome,i); //añade el sketch inicializado con 0s

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
	(kmer>comp) ? hll->insert(comp,i) : hll->insert(kmer,i); //insert kmer canonico

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

				(kmer>comp) ? hll->insert(comp,i) : hll->insert(kmer,i); //insert kmer canonico
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
	//printf("archivo: %s\n",paths);
	ifstream indata(paths);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",paths);
		exit(1);
	}
	vector<string> genomes;
	string filename;
	while(!indata.eof()){
		getline(indata,filename);
		genomes.push_back(filename);
	}
	if(genomes.back()=="") genomes.pop_back(); //elimina el EOF
	indata.close();
	return genomes;
}

//obtiene archivos de la linea de argumentos
vector<string> getPaths(char** argv, int argc){
	vector<string> genomes;
	for(int i=1;i<argc;++i){
		if(!strcmp(argv[i],"-k") || !strcmp(argv[i],"-p")) ++i;
		else if(!strcmp(argv[i],"-s")) genomes.push_back(argv[i]);
	}
	return genomes;
}

//formato: ./hll -opcion valor genomas
//o bien ./hll genomas -opcion valor
//no detecta caso en que se introduza opcion o valor invalido
int main(int argc, char *argv[]){
	if(argc<3) {
		printf("No hay suficientes argumentos");
		exit(1);
	}
	unsigned char p=22;
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
	vector<string> genomes;
	option=std::find((char**)argv,end,(const std::string&)"-f"); //lee de txt
	if(option!=end) genomes=readFromFile(*(option+1));
	else genomes=getPaths(argv,argc);
	int tam=genomes.size();

	for(vector<string>::iterator it=genomes.begin();it!=genomes.end();++it)
		printf("%s ",(char*)(*it).c_str());

	printf("\nk: %d p: %d\n",k,p);
	

	int numThreads=min(tam,(int)std::thread::hardware_concurrency());

	option=std::find((char**)argv,end,(const std::string&)"-t"); //guarda los sketches
	if(option!=end) numThreads=atoi(*(option+1));

	printf("threads: %d\n",numThreads);
	hll = new HyperLogLog(p,32-p,k,numThreads,tam);

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
				leer((char*)genomes[i].c_str(),i);
			}
		}
	}
	option=std::find((char**)argv,end,(const std::string&)"-s"); //guarda los sketches
	if(option!=end) hll->saveSketches();
	hll->estCard();
	hll->estJaccard();
	delete hll;
	return 0;
}
