BINARIES := hll
all: $(BINARIES)

CFLAGS := $(CFLAGS) -lm -pthread -Ofast -fopt-info-vec -funroll-loops -frename-registers -fno-signed-zeros -fno-trapping-math -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops  

clean:
	rm -f *.o $(BINARIES)

tags:
	etags *.h *.c *.cc

%.o: %.cpp
	g++ -c $(CFLAGS) $< -o $@ -lm -pthread -Ofast -fopt-info-vec -funroll-loops -frename-registers -fno-signed-zeros -fno-trapping-math -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops  

%.o: %.cpp
	gcc -c $(CFLAGS) $< -o $@ -lm -pthread -Ofast -fopt-info-vec -funroll-loops -frename-registers -fno-signed-zeros -fno-trapping-math -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops  

hll: main.o HyperLogLog.o  
	g++ $(CFLAGS) $^ -o $@ -lm -pthread -Ofast -fopt-info-vec -funroll-loops -frename-registers -fno-signed-zeros -fno-trapping-math -march=native -flto -fopenmp -D_GLIBCXX_PARALLEL -fipa-bit-cp -ftree-loop-vectorize -ftree-slp-vectorize -ftracer -fsplit-loops -funswitch-loops  
