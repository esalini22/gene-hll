#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H
#include <vector>
#include <cstring>
typedef unsigned long long int ullint;
using namespace std;
class HyperLogLog{
    private:
        unsigned int p;
        unsigned int b;
        unsigned int bits_v1,bits_v2; //bits para hacer & en operacion bitwise
        long double N; //corresponde a |M|: numero de buckets - se usa long double para evitar hacer casting al calcular cardinalidad
        unsigned int cerosA,cerosB,cerosU;
        vector<unsigned char> sketchA,sketchB; //se usa char pues numero más grande que puede almacenenar es 31, que usa 5 bits (1 byte = 8 bits)
        vector<unsigned int> wA,wB,wU; //vectores de contadores de w (elementos en sketches)
        pair<ullint,ullint> fractionA,fractionB,fractionU; //pares numerador-denominador
    public:
        HyperLogLog();
        ~HyperLogLog();
        void insertA(ullint kmer);
        void insertB(ullint kmer);
        void estJaccard();
};
#endif