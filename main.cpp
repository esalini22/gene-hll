#include <stdio.h>
#include <fstream>
#include "HyperLogLog.h"
using namespace std;
typedef unsigned long long int ullint;

//funcion que imprime los kmers como strings con sus respectivas bases
/*void to_kmer(unsigned long long int num,uint_fast8_t k){
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
	cout<<"kmer:"<<kmer<<endl;
}*/

int main(int argc, char *argv[]){
	//se trabaja a nivel de bits, es mas rapido que trabajar con strings
	HyperLogLog *hll = new HyperLogLog(21,11);
	char c;
	unsigned char k=31; //20 y 30
	const ullint desp=(2*(k-1));
	const ullint bits_C=(ullint)1<<desp;
	const ullint BITS=(bits_C-1)<<2;
	const ullint bits_G=(ullint)2<<desp;
	const ullint bits_T=(ullint)3<<desp;
	ullint kmer=0,comp=0;
	ullint cardA,cardB,cardU;

	//se lee linea por linea
	//se la linea en memoria para ser leida caracter por caracter

	//lee primer archivo
	string linea;
	ifstream indata(argv[1]);
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
	kmer=comp=0;

	indata.open(argv[2]);
	getline(indata,linea); //salta primera linea
	getline(indata,linea); //para primer kmer
	it=linea.begin();

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
	
	hll->estJaccard();
	delete hll;
	return 0;
}
