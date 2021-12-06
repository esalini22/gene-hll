#include "HyperLogLog.h"
#include <cmath>
#include <bitset>
#include <cstdio>
#include "wyhash32.h"

using namespace std;
#define seed 5 //seed para hash
#define lim 4294967296 //2^32

HyperLogLog::HyperLogLog(unsigned char n1, unsigned char n2){
	p=n1;
	b=n2;
	for(unsigned char i=0;i<b+1;++i){
		wA.emplace_back(0);
		wB.emplace_back(0);
		wU.emplace_back(0);
	}
	bits_v2=(1<<b)-1; //2^V2 -1
	N=1<<p; //2^V1
	a_m=(0.7213/(1+(1.079/N)))*N*N;
	cerosA=cerosB=cerosU=N;
	for(unsigned int i=0;i<N;++i){
		sketchA.emplace_back(0);
		sketchB.emplace_back(0);
	}
}
HyperLogLog::~HyperLogLog(){}

void HyperLogLog::insertA(ullint kmer){
	//calculamos el hash usando wyhash de 32 bits (hashing muy rapido)
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	v1=hash>>b; //se desplaza el hash V2 bits a la derecha
	v2=hash&bits_v2; // AND 2^V2 -1
	v2=__builtin_clz(v2)+1-p; //cuenta cantidad de bits 0 mas significativo
	//lo cual sirve para encontrar posicion de 1 mas a la izquierda
	//se resta V1, ya que al trabajar con int, se tiene 32 bits,
	//y siempre habran al menos V1 ceros a la izquierda
	if(sketchA[v1]==0) --cerosA; //bucket ya no esta en 0
	if(v2 > sketchA[v1]) sketchA[v1]=v2; //se queda con mayor valor
}
void HyperLogLog::insertB(ullint kmer){
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	v1=hash>>b;
	v2=hash&bits_v2; // AND 2^V2 -1
	v2=__builtin_clz(v2)+1-p; //posicion de primer 1, de izq a der
	if(sketchB[v1]==0) --cerosB; //bucket ya no esta en 0
	if(v2 > sketchB[v1]) sketchB[v1]=v2; //se queda con mayor valor
}
void HyperLogLog::estJaccard(){
	vector<unsigned char>::iterator it2=sketchB.begin();
	//contamos las repeticiones de w en cada sketch
	for(vector<unsigned char>::iterator it=sketchA.begin();it!=sketchA.end();++it){ 
		//se trabaja con fracciones, pues los valores son demasiado pequeños para almacenar en flotante
		unsigned char i1=*it,i2=*it2;
		wA[i1]++;
		wB[i2]++;
		if(i1||i2) --cerosU;
		(i1>i2) ? wU[i1]++ : wU[i2]++;

		++it2; //avanza en tabla B
	}
	//la idea dejar el numerador como la potencia mas grande
	//y sumar los denominadores como el contador 
	//multiplicado por la potencia de 2 de la diferencia entre el denominador original y el mayor
	unsigned char maxA=0,maxB=0,maxU=0;
	for(char i=b;i>=0;i--){ //encuentra el valor de w mas grande con cont>0
		if(maxA==0 && wA[i]) maxA=i;
		if(maxB==0 && wB[i]) maxB=i;
		if(maxU==0 && wU[i]) maxU=i;
		if(maxA && maxB && maxU) break;
	}
	fractionA.second=(ullint)1<<maxA; //2^maxA
	fractionB.second=(ullint)1<<maxB; //2^maxB
	fractionU.second=(ullint)1<<maxU; //2^maxU
	fractionA.first=fractionB.first=fractionU.first=0;
	for(unsigned char i=0;i<b;++i){
		if(wA[i]) fractionA.first+=wA[i]*((ullint)1<<(maxA-i));
		if(wB[i]) fractionB.first+=wB[i]*((ullint)1<<(maxB-i));
		if(wU[i]) fractionU.first+=wU[i]*((ullint)1<<(maxU-i));
	}
	
	//estimamos las cardinalidades
	long double cardA,cardB,cardU;
	cardA=a_m*fractionA.second/fractionA.first; //media armonica
	cardB=a_m*fractionB.second/fractionB.first;
	cardU=a_m*fractionU.second/fractionU.first;
	if(cerosA && cardA<=5*N/2) cardA=N*log(N/cerosA); //C_HLL, ln cuando hay muchos ceros;
	else if(cardA>lim/30) cardA=-lim*log(1-(cardA/lim)); //cuando la estimacion es muy grande
	printf("estimacion cardinalidad A: %Lf\n",cardA);

	if(cerosB && cardB<=5*N/2) cardB=N*log(N/cerosB);
	else if(cardB>lim/30) cardB=-lim*log(1-(cardB/lim));
	printf("estimacion cardinalidad B: %Lf\n",cardB);

	if(cerosU && cardU<=5*N/2) cardU=N*log(N/cerosU);
	else if(cardU>lim/30) cardU=-lim*log(1-(cardU/lim));
	printf("estimacion cardinalidad A U B: %Lf\n",cardU);

	printf("estimacion jaccard: %Lf\n",(cardA+cardB-cardU)/cardU);
}
