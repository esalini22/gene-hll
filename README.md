# HyperLogLog-cpp
HyperLogLog en c++ para cálculo de distancia de genomas mediante índice de Jaccard.

Hace uso de la función de hashing wyhash32, pues se encontró que es la más rápida, sin afectar negativamente a los resultados.

El algoritmo trata a los k-mers como cadenas de bits, donde cada base es representada por 2 bits, con A=00, C=01, G=10, T=11. Como el largo máximo del kmer es de 31, utiliza a lo más 62 bits, lo cual cabe en un long long.


## Uso
Se compila de la siguiente forma:

g++ -c main.cpp HyperLogLog.cpp -Ofast -fopt-info-vec -funroll-loops -frename-registers -fno-signed-zeros -fno-trapping-math -march=native -flto -fprofile-generate

g++ -o output_file main.o HyperLogLog.o -Ofast -fopt-info-vec -funroll-loops -frename-registers -fno-signed-zeros -fno-trapping-math -march=native -flto -fprofile-generate

Se ejecuta una vez, y se vuelve a compilar, reemplazando -fprofile-generate por -fprofile-use

Se ejecuta como *./output_file genoma1 genoma2*

## Comparación con Dashing
