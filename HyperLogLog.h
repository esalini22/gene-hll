#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H
#include <vector>
//#include <unordered_map>
#include <string>
typedef unsigned long long int ullint;
using namespace std;
class HyperLogLog{
    private:
        unsigned char p,b;
        unsigned int bits_v2; //bits para hacer & en operacion bitwise
        long double N; //corresponde a |M|: numero de buckets - se usa long double para evitar hacer casting al calcular cardinalidad
        long double a_m; //factor de correccion * N^2
        string kmer_length;
        string sketch_size;
        vector<ullint> bit_mask; //vector de máscaras de bits, usadas para insertar registros
        //se usa vector de ullint, cada celda (64 bits) a su vez tiene hasta 12 buckets
        vector<vector<ullint>> sketch; //n sketches, n: numero de genomas
        vector<string> names; //vector con los nombres de los genomas
        vector<long double> cards; //cardinalidades individuales, para cálculo de jaccard
        //unordered_map<string,long double> jaccards; //mapa con jaccard de pares, para matriz de dist
    public:
        HyperLogLog(unsigned char n1,unsigned char n2,unsigned char k,int numThreads,int n);
        ~HyperLogLog();
        void addSketch(string genome,int i); //añade un sketch vacio al mapa de sketches, para funciones hll y dist
        void insert(ullint kmer,int i); //inserta en sketch, para funciones hll y dist
        void loadSketch(char* filename); //carga un archivo de sketch en memoria, para hll y dist
        void saveSketches(); //guarda el sketch individual en un archivo, para funcion sketch
        void estCard(); //estima cardinalidad individual, para funciones hll y dist
        void estJaccard(); //estima jaccard, para funcion dist
};
#endif
