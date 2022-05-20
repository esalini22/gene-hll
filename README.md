# gene-hll
Implementación del algoritmo HyperLogLog en C++ para cálculo de similitud de genomas mediante índice de Jaccard.
El algoritmo se basa en el siguiente [paper](https://storage.googleapis.com/pub-tools-public-publication-data/pdf/40671.pdf).

Recibe como entrada dos archivos de genomas, y devuelve:
- Las cardinalidades de cada genoma
- El índice de Jaccard de cada par de genomas

Hace uso de la función de hashing [wyhash32](https://github.com/wangyi-fudan/wyhash), pues se encontró que es la más rápida, sin afectar negativamente a los resultados.

## Uso
Se compila usando make

Se ejecuta como *./hll -opcion valor genomas*, o bien *./hll genomas -opcion valor*

### Opciones
- -p: Cambia el valor de los p bits. Mínimo: 7, máximo: 31
- -k: Cambia el largo del kmer. Mínimo: 20, máximo: 31
- -f: Lee las rutas de los archivos de genomas desde un archivo de texto.
- -s: Guarda los sketches generados en un archivo de texto.
- -t: Número de hebras

Ejemplo:

*./hll -p 22 genoma1 genoma 2*
*./hll -k 31 -f archivo*
*./hll -s genoma1 genoma 2 genoma3 -k 32

Nota: El código no detecta opciones inválidas.
