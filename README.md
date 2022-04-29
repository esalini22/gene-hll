# gene-hll
Implementación del algoritmo HyperLogLog en C++ para cálculo de similitud de genomas mediante índice de Jaccard.
El algoritmo se basa en el siguiente [paper](https://storage.googleapis.com/pub-tools-public-publication-data/pdf/40671.pdf).

Recibe como entrada dos archivos de genomas, y devuelve:
- Las cardinalidades de cada genoma
- La cardinalidad de la unión
- El índice de Jaccard

Hace uso de la función de hashing [wyhash32](https://github.com/wangyi-fudan/wyhash), pues se encontró que es la más rápida, sin afectar negativamente a los resultados.

## Uso
Se compila usando make

Se ejecuta como *./hll genoma1 genoma2*

### Opciones
- p: Cambia el valor de los p bits. Mínimo: 7, máximo: 31
- k: Cambia el largo del kmer. Mínimo: 20, máximo: 31

Ejemplo:

*./hll genoma1 genoma2 k 20 p 14*
