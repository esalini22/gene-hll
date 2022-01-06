BINARIES := hll
all: $(BINARIES)

CFLAGS := $(CFLAGS) -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 

clean:
	rm -f *.o $(BINARIES)

tags:
	etags *.h *.c *.cc

%.o: %.cpp
	g++ -c $(CFLAGS) $< -o $@ -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 

%.o: %.cpp
	gcc -c $(CFLAGS) $< -o $@ -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 

hll: main.o HyperLogLog.o  
	g++ $(CFLAGS) $^ -o $@ -lm -pthread -Ofast -funroll-loops -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 
