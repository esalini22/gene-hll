#include "HyperLogLog.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include <zlib.h>
#include <stdio.h>
#include "wyhash32.h"

using namespace std;
#define seed 5 //seed para hash
#define lim 4294967296 //2^32
#define lim_int 2147483647

HyperLogLog::HyperLogLog(unsigned char n1, unsigned char n2,unsigned char k,int numThreads, int n){
	omp_set_num_threads(numThreads);
	p=n1;
	b=n2;
	bits_v2=(1<<b)-1; //2^b -1
	N=1<<p; //2^V1
	a_m=(0.7213/(1+(1.079/N)))*N*N;
	for(unsigned char i=0;i<16;++i)
		bit_mask.emplace_back(~((ullint)0xF<<(4*i)));
	cards.resize(n);
	sketch.resize(n);
	names.resize(n);
	kmer_length=to_string((int)k);
	sketch_size=to_string(p);
	jaccards.resize(n);
	ciclos_red=(b+2)/8+(((b+2)%8)>0);
	min_val=1<<(b-15);
}
HyperLogLog::~HyperLogLog(){}

void HyperLogLog::addSketch(string genome,int i){
	//inicializa sketch en 0s, y añade su nombre a un vector
	for(unsigned int j=0;j<N/16;++j) //cada celda tendra 12 buckets del sketch array
		sketch[i].emplace_back(0);
	if((long)N%16)
		sketch[i].emplace_back(0);
	names[i]=genome;
}

void HyperLogLog::insert(ullint kmer,int i){
	//calculamos el hash usando wyhash de 32 bits (hashing muy rapido)
	const ullint *key = &kmer;
	uint32_t hash=wyhash32(key,8,seed);

	unsigned int v1=hash>>b; //se desplaza el hash b bits a la derecha
	unsigned int v2=hash&bits_v2; // AND 2^b -1 
	if(v2<min_val) v2=15;
	else v2=__builtin_clz(v2)+1-p; //cuenta cantidad de bits 0 mas significativo
	//lo cual sirve para encontrar posicion de 1 mas a la izquierda
	//se resta p, ya que al trabajar con int, se tiene 32 bits,
	//y siempre habran al menos p ceros a la izquierda
	unsigned int indice=v1/16; //indice de celda
	unsigned char bucket=v1%16; //12 buckets por celda
	ullint temp=(sketch[i][indice]>>(4*bucket))&0xF;
	if(v2 > temp) sketch[i][indice]=(sketch[i][indice]&(bit_mask[bucket]))|((ullint)v2<<(4*bucket));
}

void HyperLogLog::loadSketch(char *infilename,int i){ //¿que pasa si tamaño del sketch es distinto al de argumento?
	string genome=infilename; //se debe determinar nombre del genoma a partir del archivo .hll
	printf("infilename: %s\n",infilename);
	genome.erase(genome.find(".w."),string::npos);
	printf("genome: %s\n",genome.c_str());

	//debemos descomprimir el archivo y ponerlo en un vector
	vector<ullint> v;
	gzFile infile = gzopen(infilename, "rb");
	//if(!infile) return -1;
	char buffer[128];
	int num_read = 0;
	char temp[19];
	char* end;
	int cont=0;
	while((num_read = gzread(infile,buffer, sizeof(buffer))) > 0){
		for(int i=0;i<num_read;++i){
			temp[cont%19]=buffer[i];
			++cont;
			if(cont%19==0){
				cont=0;
				v.emplace_back(strtoull(temp,&end,10));
			}
		}
	}
	if(cont>0) v.emplace_back(strtoull(temp,&end,10));
	v.emplace_back(0);
	gzclose(infile);

	sketch[i]=v;
	names[i]=genome;
	printf("done\n");
}

void HyperLogLog::saveSketches(){ //usar zlib para comprimir
	printf("save sketches\n");
	int tam=sketch.size();
	#pragma omp parallel for
	for(int i=0;i<tam;++i){
		vector<ullint> s=sketch[i];
		string name=names[i];
		name+=".w."+kmer_length+".spacing."+sketch_size+".hll";
		char* outfilename=(char*)name.c_str();

		printf("%s\n",outfilename);

		gzFile outfile = gzopen(outfilename, "wb");
		//if(s.empty() || !outfile) return -1;
		char inbuffer[128];
		char cont=0;
		for(vector<ullint>::iterator it=s.begin();it!=s.end();++it){
			string temp=to_string(*it)+'\n';
			int len=temp.length();
			for(char i=0;i<len;++i){
				inbuffer[cont]=temp[i];
				++cont;
				if(cont==127){
					gzwrite(outfile, inbuffer, cont);
					cont=0;
				}
			}
		}
		if(cont>0) gzwrite(outfile, inbuffer, cont);
		gzclose(outfile);
	}
}

void HyperLogLog::saveOutput(char* filename){ //guarda la matriz en txt
	int tam=sketch.size();
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

//siguientes dos funciones son para reduce vectorizado

inline float HyperLogLog::hsum_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 maxs = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, maxs); // high half -> low half
    maxs        = _mm_add_ss(maxs, shuf);
    return        _mm_cvtss_f32(maxs);
}

inline float HyperLogLog::hsum_avx(__m256 v) {
    __m128 lo = _mm256_castps256_ps128(v);   // low 128
    __m128 hi = _mm256_extractf128_ps(v, 1); // high 128
           lo = _mm_add_ps(lo, hi);          // max the low 128
    return hsum_sse3(lo);                    // and inline the sse3 version
}

void HyperLogLog::estCard(){
	__m256i vec3; //vector de mascaras de bits
	ullint bits_and[4]={15,15,15,15};
	//vec3=_mm256_loadu_epi64((void const *)&bits_and[0]);
	vec3=_mm256_loadu_si256((const __m256i*)&bits_and[0]);

	int tam=sketch.size();
	#pragma omp parallel for
	for (int j = 0; j < tam; ++j){
		FILE *fp = fopen((names[j]+".sketch").c_str(),"w");

		int ceros=N;
		vector<float> w(b+2,0);
		vector<ullint>::iterator it1=sketch[j].begin();
		vector<ullint>::iterator fin=sketch[j].end();
		//contamos las repeticiones de w en cada sketch
		while(it1!=fin){ 
			ullint i1=*it1;

			ullint it_array[16]={i1,(i1>>4),(i1>>8),(i1>>12),(i1>>16),(i1>>20),(i1>>24),(i1>>28),(i1>>32),(i1>>36),(i1>>40),(i1>>44),(i1>>48),(i1>>52),(i1>>56),(i1>>60)};
			__m256i vec4[4];
			vec4[0]=_mm256_loadu_si256((const __m256i*)&it_array[0]);
			vec4[1]=_mm256_loadu_si256((const __m256i*)&it_array[4]);
			vec4[2]=_mm256_loadu_si256((const __m256i*)&it_array[8]);
			vec4[3]=_mm256_loadu_si256((const __m256i*)&it_array[12]);
			vec4[0]=_mm256_and_si256(vec3,vec4[0]);
			vec4[1]=_mm256_and_si256(vec3,vec4[1]);
			vec4[2]=_mm256_and_si256(vec3,vec4[2]);
			vec4[3]=_mm256_and_si256(vec3,vec4[3]);
			__attribute__ ((aligned (32))) ullint out[4];
			for(int c=0;c<4;++c){
				_mm256_store_si256((__m256i *)&out[0],vec4[c]);
				if(out[0]) ceros--;
				w[out[0]]++;
				if(out[1]) ceros--;
				w[out[1]]++;
				if(out[2]) ceros--;
				w[out[2]]++;
				if(out[3]) ceros--;
				w[out[3]]++;
			}
			
			++it1;
		}

		float card=0;
		vector<float> w2;
		for(int i=0;i<b+2;++i) w2.emplace_back((float)1/(float)(1<<i));
		__m256 vec,vec2;
		for(int i=0;i<ciclos_red;++i){
			vec=_mm256_loadu_ps((const float *)&w[i*8]);
			vec2=_mm256_loadu_ps((const float *)&w2[i*8]);
			vec=_mm256_mul_ps(vec,vec2);
			card+=hsum_avx(vec);
		}

		//media armonica
		card=(float)a_m/card;
		if(ceros && card<=5*N/2) //C_HLL, ln cuando hay muchos ceros;
			card=N*log(N/ceros);
		else if(card>lim/30)
			card=-lim*log(1-(card/lim));
		printf("estimacion cardinalidad %s: %f\n",names[j].c_str(),card);
		cards[j]=card;

		fclose(fp);
	}
	printf("\n");
}

void HyperLogLog::estJaccard(){ //recibe numero de genomas
	__m256i vec3; //vector de mascaras de bits
	ullint bits_and[4]={15,15,15,15};
	vec3=_mm256_loadu_si256((const __m256i*)&bits_and[0]);

	int tam=sketch.size();
	#pragma omp parallel for //no se puede usar collapse aqui
	for(int j1=0;j1<tam;j1++){
		vector<float> temp_jaccards; //para evitar false sharing
		temp_jaccards.resize(tam-j1-1);
		for(int j2=j1+1;j2<tam;j2++){
			vector<ullint>::iterator it1=sketch[j1].begin();
			vector<ullint>::iterator it2=sketch[j2].begin();
			vector<ullint>::iterator fin=sketch[j1].end();
			int ceros=N;
			vector<float> wU(b+2,0);
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
					if(temp1||temp2) --ceros;
					(temp1>temp2) ? wU[temp1]++ : wU[temp2]++;
					i1=i1>>4;
					i2=i2>>4;
				}
				++it1;
				++it2; //avanza en tabla B
			}
			
			float cardU=0;

			vector<float> w2;
			for(int i=0;i<min(b+2,17);++i) w2.emplace_back((float)1/(float)(1<<i));
			__m256 vec,vec2;
			for(int i=0;i<ciclos_red;++i){
				vec=_mm256_loadu_ps((const float *)&wU[i*8]);
				vec2=_mm256_loadu_ps((const float *)&w2[i*8]);
				vec=_mm256_mul_ps(vec,vec2);
				cardU+=hsum_avx(vec);
			}

			//media armonica
			cardU=(float)a_m/cardU;
			if(ceros && cardU<=5*N/2) //C_HLL, ln cuando hay muchos ceros;
				cardU=N*log(N/ceros);
			else if(cardU>lim/30)
				cardU=-lim*log(1-(cardU/lim));
			//printf("estimacion cardinalidad %s U %s: %Lf\n",names[j1].c_str(),names[j2].c_str(),cardU);

			float jaccard=(cards[j1]+cards[j2]-cardU)/cardU;
			if(jaccard<0) jaccard=0;
			temp_jaccards[j2-j1-1]=jaccard;
		}
		jaccards[j1]=temp_jaccards;
	}
}

void HyperLogLog::printMatrix(){
	int tam=sketch.size();
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
