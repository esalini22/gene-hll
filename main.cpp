#include <iostream>
#include <sys/mman.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <cmath>
#include "HyperLogLog.h"
using namespace std;
typedef unsigned long long int ullint;

/*#define BITS 0xFFFFFFFFFC
#define bits_C 0x4000000000
#define bits_G 0x8000000000
#define bits_T 0xC000000000*/

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
	const char *file1 = argv[1], *file2 = argv[2];
	HyperLogLog *hll = new HyperLogLog();
	char c;
	unsigned char k=31; //20 y 30
	const ullint desp=(2*(k-1));
	const ullint bits_C=(ullint)1<<desp;
	const ullint BITS=(bits_C-1)<<2;
	const ullint bits_G=(ullint)2<<desp;
	const ullint bits_T=(ullint)3<<desp;
	ullint kmer=0,comp=0;

	//lee mapeando archivo a memoria
	//lee primer archivo
	int fd = open(file1, O_RDONLY);
	struct stat stat;
	int err = fstat(fd, &stat);
    char *ptr = (char*)mmap(NULL,stat.st_size,PROT_READ,MAP_PRIVATE | MAP_POPULATE,fd,0);
    uint_fast32_t pos=0;

	while((c=*(ptr+pos))!='\n') ++pos; //lee primera linea

	for(unsigned char i=0;i<k;++i){ //inicializacion - primer kmer
		++pos;
		kmer=kmer<<2;
		comp=comp>>2;
		c = *(ptr+pos);
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
	while(c!=EOF){ //genera kmer y lo inserta a hll
		++pos;
		c = *(ptr+pos);
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
		else if(c=='>') while((c=*(ptr+pos))!='\n') ++pos; //se salta la linea
	}
	close(fd);
	err = munmap(ptr, stat.st_size);

	//lee segundo archivo
	fd = open(file2, O_RDONLY);
	err = fstat(fd, &stat);
	ptr = (char*)mmap(NULL,stat.st_size,PROT_READ,MAP_PRIVATE | MAP_POPULATE,fd,0);
	pos=0;
	kmer=0;
	comp=0;

	while((c=*(ptr+pos))!='\n') ++pos; //lee la primera linea

	for(unsigned char i=0;i<k;++i){ //inicializacion - primer kmer
		++pos;
		kmer=kmer<<2;
		comp=comp>>2;
		c = *(ptr+pos);
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
	(kmer>comp) ? hll->insertB(comp) : hll->insertB(kmer); //inserta kmer canonico
	while(c!=EOF){ //genera kmer y lo inserta a hll
		++pos;
		c = *(ptr+pos);
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
			(kmer>comp) ? hll->insertB(comp) : hll->insertB(kmer); //inserta kmer canonico
		}
		else if(c=='>') while((c=*(ptr+pos))!='\n') ++pos; //se salta la linea
	}
	close(fd);
	err = munmap(ptr, stat.st_size);
	
	hll->estJaccard();
	delete hll;
	return 0;
}