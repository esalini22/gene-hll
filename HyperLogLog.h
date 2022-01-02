#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H
#include <vector>
typedef unsigned long long int ullint;
using namespace std;
class HyperLogLog{
    private:
        unsigned char p,b;
        unsigned int bits_v2; //bits para hacer & en operacion bitwise
        long double N; //corresponde a |M|: numero de buckets - se usa long double para evitar hacer casting al calcular cardinalidad
        long double a_m; //factor de correccion * N^2
        unsigned int cerosA,cerosB,cerosU;
        vector<ullint> bit_mask; //vector de máscaras de bits, usadas para insertar registros
        vector<ullint> sketchA,sketchB; //se usa vector de ullint, cada celda (64 bits) a su vez tiene hasta 12 buckets
        vector<unsigned int> wA,wB,wU; //vectores de contadores de w (elementos en sketches)
    public:
        HyperLogLog(unsigned char n1,unsigned char n2);
        ~HyperLogLog();
        void insertA(ullint kmer);
        void insertB(ullint kmer);
        void estJaccard();
};
#endif
