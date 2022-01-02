#include "HyperLogLog.h"
#include <cmath>
#include <bitset>
#include <cstdio>
#include "wyhash32.h"

using namespace std;
#define seed 5
#define lim 4294967296 //2^32

HyperLogLog::HyperLogLog(unsigned char n1, unsigned char n2){
	p=n1;
	b=n2;
	for(unsigned char i=0;i<b+2;++i){ //es b+2 en vez de b+1, pues v2 podria ser 0 (PREGUNTAR)
		wA.emplace_back(0);
		wB.emplace_back(0);
		wU.emplace_back(0);
	}
	bits_v2=(1<<b)-1; //2^V2 -1
	N=1<<p; //2^V1
	a_m=(0.7213/(1+(1.079/N)))*N*N;
	cerosA=cerosB=cerosU=N;
	for(unsigned int i=0;i<N/12;++i){ //cada celda tendra 12 buckets del sketch array
		sketchA.emplace_back(0);
		sketchB.emplace_back(0);
	}
	if((long)N%12){
		sketchA.emplace_back(0);
		sketchB.emplace_back(0);
	}
}
HyperLogLog::~HyperLogLog(){}

void HyperLogLog::insertA(ullint kmer){
	//calculamos el hash usando wyhash de 32 bits (hashing muy rapido)
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	unsigned int v1=hash>>b; //se desplaza el hash V2 bits a la derecha
	ullint v2=hash&bits_v2; // AND 2^V2 -1
	if(v2==0) v2=b+1;
	else v2=__builtin_clz(v2)+1-p; //cuenta cantidad de bits 0 mas significativo
	//lo cual sirve para encontrar posicion de 1 mas a la izquierda
	//se resta p, ya que al trabajar con int, se tiene 32 bits,
	//y siempre habran al menos p ceros a la izquierda
	unsigned int indice=v1/12; //indice de celda
	unsigned char bucket=v1%12; //12 buckets por celda
	ullint temp=(sketchA[indice]>>(5*bucket))&0x1F;
	if(v2 > temp){
		if(temp==0) --cerosA; //bucket ya no esta en 0
		ullint bit_mask=(ullint)0x1F<<(5*bucket);
		sketchA[indice]=(sketchA[indice]&(~bit_mask))|(v2<<(5*bucket));
	}
}
void HyperLogLog::insertB(ullint kmer){
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	unsigned int v1=hash>>b;
	ullint v2=hash&bits_v2; // AND 2^V2 -1
	if(v2==0) v2=b+1;
	else v2=__builtin_clz(v2)+1-p; //posicion de primer 1, de izq a der

	unsigned int indice=v1/12;
	unsigned char bucket=v1%12;
	ullint temp=(sketchB[indice]>>(5*bucket))&0x1F;

	if(v2 > temp){ //se queda con mayor valor
		if(temp==0) --cerosB; //bucket ya no esta en 0
		ullint bit_mask=(ullint)0x1F<<(5*bucket);
		sketchB[indice]=(sketchB[indice]&(~bit_mask))|(v2<<(5*bucket));
	}
}
void HyperLogLog::estJaccard(){
	long double cardA=0,cardB=0,cardU=0;
	vector<ullint>::iterator it1=sketchA.begin();
	vector<ullint>::iterator it2=sketchB.begin();
	vector<ullint>::iterator fin=sketchA.end();
	while(it1!=fin){
		ullint i1=*it1,i2=*it2;
		for(char i=0;i<12;++i){ //12 registros por celda
			ullint temp1=i1&0x1F,temp2=i2&0x1F;
			wA[temp1]++;
			wB[temp2]++;
			if(temp1||temp2) --cerosU;
			(temp1>temp2) ? wU[temp1]++ : wU[temp2]++;
			i1=i1>>5;
			i2=i2>>5;
		}
		++it1;
		++it2; //avanza en tabla B
	}
	for(unsigned char i=0;i<b+2;++i){
		if(wA[i]) cardA+=(long double)wA[i]/(long double)(1<<i);
		if(wB[i]) cardB+=(long double)wB[i]/(long double)(1<<i);
		if(wU[i]) cardU+=(long double)wU[i]/(long double)(1<<i);
	}
	//media armonica
	cardA=(long double)a_m/cardA;
	cardB=(long double)a_m/cardA;
	cardU=(long double)a_m/cardA;
	if(cerosA && cardA<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
		//printf("ceros A:%u\n",cerosA);
		printf("linear counting\n");
		cardA=N*log(N/cerosA);
	}
	else if(cardA>lim/30){
		printf("valor muy grande\n");
		cardA=-lim*log(1-(cardA/lim));
	}
	printf("estimacion cardinalidad A: %Lf\n",cardA);

	if(cerosB && cardB<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
		//printf("ceros B:%u\n",cerosB);
		printf("linear counting\n");
		cardB=N*log(N/cerosB);
	}
	else if(cardB>lim/30){
		printf("valor muy grande\n");
		cardB=-lim*log(1-(cardB/lim));
	}
	printf("estimacion cardinalidad B: %Lf\n",cardB);

	if(cerosU && cardU<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
		//printf("ceros A U B:%u\n",cerosU);
		printf("linear counting\n");
		cardU=N*log(N/cerosU);
	}
	else if(cardU>lim/30){
		printf("valor muy grande\n");
		cardU=-lim*log(1-(cardU/lim));
	}
	printf("estimacion cardinalidad A U B: %Lf\n",cardU);

	printf("estimacion jaccard: %Lf\n",(cardA+cardB-cardU)/cardU);
}
