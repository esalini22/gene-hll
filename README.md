# HyperLogLog-cpp
HyperLogLog en c++ para cálculo de distancia de genomas mediante índice de Jaccard.

Hace uso de la función de hashing wyhash32, pero puede ser reemplazada por cualquier otra de 32bits.

Trata a los k-mers como cadenas de bits, donde cada base es representada por 2 bits. Donde A=00, C=01, G=10, T=11.
