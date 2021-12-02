#include "HyperLogLog.h"
#include <cmath>
#include <bitset>
#include <cstdio>
#include "wyhash32.h"

using namespace std;
//V1+V2=32
#define V1 14 //determinar en base a error = 1.04/sqrt(|M|)
#define V2 18
#define seed 5 //seed para hash - CAMBIAR PARA VER COMO CAMBIA TIEMPO Y PRECISION
#define lim 4294967296 //2^32
//#define N 512
//#define bits_v2 8388607

HyperLogLog::HyperLogLog(){
	N=1<<V1; //2^V1
	cerosA=N;
	cerosB=N;
	cerosU=N;
	bits_v2=(1<<V2)-1; //2^V2 -1
	sketchA.clear();
	sketchB.clear();
	sketchA.resize(N);
	sketchB.resize(N);
	wA.clear();
	wB.clear();
	wU.clear();
	wA.resize(V2+1);
	wB.resize(V2+1);
	wU.resize(V2+1);
}
HyperLogLog::~HyperLogLog(){}

void HyperLogLog::insertA(ullint kmer){
	//calculamos el hash usando wyhash de 32 bits (hashing muy rapido)
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	p=hash>>V2; //se desplaza el hash V2 bits a la derecha
	b=hash&bits_v2; // AND 2^V2 -1
	//b=(hash<<V1)>>V1;
	b=__builtin_clz(b)+1-V1; //cuenta cantidad de bits 0 mas significativo
	//lo cual sirve para encontrar posicion de 1 mas a la izquierda
	//se resta V1, ya que al trabajar con int, se tiene 32 bits,
	//y siempre habran al menos V1 ceros a la izquierda
	if(sketchA[p]==0) --cerosA; //bucket ya no esta en 0
	if(b > sketchA[p]) sketchA[p]=b; //se queda con mayor valor
}
void HyperLogLog::insertB(ullint kmer){
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	p=hash>>V2;
	b=hash&bits_v2; // AND 2^V2 -1
	b=__builtin_clz(b)+1-V1; //posicion de primer 1, de izq a der
	if(sketchB[p]==0) --cerosB; //bucket ya no esta en 0
	if(b > sketchB[p]) sketchB[p]=b; //se queda con mayor valor
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
	for(char i=V2;i>=0;i--){ //encuentra el valor de w mas grande con cont>0
		if(maxA==0 && wA[i]) maxA=i;
		if(maxB==0 && wB[i]) maxB=i;
		if(maxU==0 && wU[i]) maxU=i;
		if(maxA && maxB && maxU) break;
	}
	fractionA.second=(ullint)1<<maxA; //2^maxA
	fractionB.second=(ullint)1<<maxB; //2^maxB
	fractionU.second=(ullint)1<<maxU; //2^maxU
	fractionA.first=0;
	fractionB.first=0;
	fractionU.first=0;
	for(unsigned char i=0;i<V2;++i){
		if(wA[i]) fractionA.first+=wA[i]*((ullint)1<<(maxA-i));
		if(wB[i]) fractionB.first+=wB[i]*((ullint)1<<(maxB-i));
		if(wU[i]) fractionU.first+=wU[i]*((ullint)1<<(maxU-i));
	}
	/*
	//simplificamos
	//solo en caso en que la multiplicacion del calculo de cardinalidades siguiente sea demasiado grande
	ullint temp1=fractionA.first,temp2=fractionA.second,mcd;
	mcd=(temp2>temp1) ? temp1:temp2;
	while(mcd>0){
		if(temp1%mcd==0 && temp2%mcd==0) break;
		--mcd;
	}
	fractionA.first=fractionA.first/mcd;
	fractionA.second=fractionA.second/mcd;

	temp1=fractionB.first,temp2=fractionB.second,mcd;
	mcd=(temp2>temp1) ? temp1:temp2;
	while(mcd>0){
		if(temp1%mcd==0 && temp2%mcd==0) break;
		--mcd;
	}
	fractionB.first=fractionB.first/mcd;
	fractionB.second=fractionB.second/mcd;

	temp1=fractionU.first,temp2=fractionU.second,mcd;
	mcd=(temp2>temp1) ? temp1:temp2;
	while(mcd>0){
		if(temp1%mcd==0 && temp2%mcd==0) break;
		--mcd;
	}
	fractionU.first=fractionU.first/mcd;
	fractionU.second=fractionU.second/mcd;*/

	/*printf("ceros A:%u\n",cerosA);
	printf("ceros B:%u\n",cerosB);
	printf("ceros U:%u\n",cerosU);*/
	long double a_m=(0.7213/(1+(1.079/N)))*N*N; //factor de correccion * N^2
	long double cardA,cardB,cardU;
	cardA=a_m*fractionA.second/fractionA.first; //media armonica
	cardB=a_m*fractionB.second/fractionB.first;
	cardU=a_m*fractionU.second/fractionU.first;
	if(cerosA && cardA<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
		printf("ceros A:%u\n",cerosA);
		printf("linear counting\n");
		cardA=N*log(N/cerosA);
	}
	else if(cardA>lim/30){
		printf("valor muy grande\n");
		cardA=-lim*log(1-(cardA/lim));
	}
	printf("estimacion cardinalidad A: %Lf\n",cardA);

	if(cerosB && cardB<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
		printf("ceros B:%u\n",cerosB);
		printf("linear counting\n");
		cardB=N*log(N/cerosB);
	}
	else if(cardB>lim/30){
		printf("valor muy grande\n");
		cardB=-lim*log(1-(cardB/lim));
	}
	printf("estimacion cardinalidad B: %Lf\n",cardB);

	if(cerosU && cardU<=5*N/2){ //C_HLL, ln cuando hay muchos ceros;
		printf("ceros A U B:%u\n",cerosU);
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