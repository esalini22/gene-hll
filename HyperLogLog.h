#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H
#include <vector>
#include <string>
#include <immintrin.h> 
typedef unsigned long long int ullint;
using namespace std;
class HyperLogLog{
    private:
        unsigned char p,b;
        unsigned int bits_v2; //bits para hacer & en operacion bitwise
        long double N; //corresponde a |M|: numero de buckets - se usa long double para evitar hacer casting al calcular cardinalidad
        float a_m; //factor de correccion * N^2
        int ciclos_red; //cantidad de ciclos de for de reduce vectorizado
        string kmer_length;
        string sketch_size;
        __m256 vect[4]; //vector wt vectorizado
        vector<ullint> bit_mask; //vector de máscaras de bits, usadas para insertar registros
        //se usa vector de ullint, cada celda (64 bits) a su vez tiene hasta 12 buckets
        vector<vector<ullint>> sketch; //n sketches, n: numero de genomas
        vector<string> names; //vector con los nombres de los genomas
        vector<float> cards; //cardinalidades individuales, para cálculo de jaccard
        vector<vector<float>> jaccards; //matriz de jaccards
        unsigned int min_val; //minimo valor en v2 para no sobrepasar registro 15
    public:
        HyperLogLog(unsigned char n1,unsigned char n2,unsigned char k,int numThreads,int n);
        ~HyperLogLog();
        void addSketch(string genome,int i); //añade un sketch vacio al mapa de sketches, para funciones hll y dist
        void insert(ullint kmer,int i); //inserta en sketch, para funciones hll y dist
        void loadSketch(char* infilename,int i); //carga un archivo de sketch en memoria, para hll y dist
        void saveSketches(); //guarda el sketch individual en un archivo, para funcion sketch
        void saveOutput(char* filename); //guarda las cardinalidades y distancias en un txt
        inline float hsum_sse3(__m128 v);
        inline float hsum_avx(__m256 v);
        void estCard(); //estima cardinalidad individual, para funciones hll y dist
        void estJaccard(); //estima jaccard, para funcion dist
        void printMatrix();
};
#endif
