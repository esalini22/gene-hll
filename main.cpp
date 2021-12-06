#include <stdio.h>
#include <fstream>
#include <unistd.h>
#include <pthread.h>
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

void* leerA(void *filename){
	char c; //para lectura
	ullint kmer=0,comp=0;
	string linea;
	ifstream indata((char*)filename);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",(char*)filename);
		exit(1);
	}
	getline(indata,linea); //salta primera linea
	getline(indata,linea); //para primer kmer
	string::iterator it=linea.begin();

	//al leerse cada caracter, se inserta la base al final del kmer no canonico
	//y luego se inserta su complemento al inicio del complemento del reverso del kmer
	//el inicio estara dado por el largo del kmer
	for(unsigned char i=0;i<k;++i){ //inicializacion - primer kmer
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
	(kmer>comp) ? hll->insertA(comp) : hll->insertA(kmer); //inserta kmer canonico

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

				(kmer>comp) ? hll->insertA(comp) : hll->insertA(kmer); //inserta kmer canonico
			}
			else if(c=='>') break;
			++it;
		}
		getline(indata,linea);
		it=linea.begin();
	}
	indata.close();
	return NULL;
}

void* leerB(void *filename){
	//lee segundo archivo
	char c; //para lectura
	ullint kmer=0,comp=0;
	string linea;
	ifstream indata((char*)filename);
	if(!indata){
		printf("No se pudo abrir el archivo %s\n",(char*)filename);
		exit(1);
	}
	getline(indata,linea); //salta primera linea
	getline(indata,linea); //para primer kmer
	string::iterator it=linea.begin();

	for(unsigned char i=0;i<k;++i){ //inicializacion - primer kmer
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
	(kmer>comp) ? hll->insertB(comp) : hll->insertB(kmer);

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
				(kmer>comp) ? hll->insertB(comp) : hll->insertB(kmer);
			}
			else if(c=='>') break;
			++it;
		}
		getline(indata,linea);
		it=linea.begin();
	}
	indata.close();
	return NULL;
}

int main(int argc, char *argv[]){
	//se trabaja a nivel de bits, es mas rapido que trabajar con strings
	if(argc<3) {
		printf("No hay suficientes argumentos");
		exit(1);
	}
	unsigned char p=21;
	if(argc==5 || argc==6){ //1 opcion: k o p
		char arg=atoi(argv[4]);
		if(!strcmp(argv[3],"k") && arg<32 && arg>19){
			printf("opcion k\n");
			k=arg;
		}
		else if(!strcmp(argv[3],"p") && arg<32 && arg>8){
			printf("opcion p\n");
			p=arg;
		}
	}
	else if(argc>6){ //2 opciones: k y p
		char arg1=atoi(argv[4]);
		char arg2=atoi(argv[6]);
		if(!strcmp(argv[3],"k") && arg1<32 && arg1>19) k=arg1;
		else if(!strcmp(argv[3],"p") && arg1<32 && arg1>8) p=arg1;
		if(!strcmp(argv[5],"k") && arg2<32 && arg2>19) k=arg2;
		else if(!strcmp(argv[5],"p") && arg2<32 && arg2>8) p=arg2;
	}
	hll = new HyperLogLog(p,32-p);
	//aca se determina como se insertaran las bases complementarias en el complemento del reverso del kmer
	//es decir, se determina donde esta el inicio (en bits) de dicho kmer
	const ullint desp=(2*(k-1));
	bits_G=(ullint)2<<desp;
	bits_T=(ullint)3<<desp;
	bits_C=(ullint)1<<desp;
	//esto sirve para eliminar la primera base del kmer, luego de desplazar las bases a la izquierda para leer la nueva base
	BITS=(bits_C-1)<<2;

	//lee paralelamente cada archivo
	pthread_t ids[2]; //arreglo de hebras
	pthread_create(&ids[0], NULL, leerA, argv[1]);
	pthread_create(&ids[1], NULL, leerB, argv[2]);
	
	pthread_join(ids[0], NULL);
	pthread_join(ids[1], NULL);
	hll->estJaccard();
	delete hll;
	return 0;
}
