CC = mpicc
MPI_FLAGS = -Wall -lm -fopenmp -O3 -g
CUDA_LIB = -L/usr/local/cuda-7.5/targets/x86_64-linux/lib -lcudart
MPI_LIB = -I/opt/tools/libraries/openmpi/2.0.1-gcc-5.4.0/include -L/opt/tools/libraries/openmpi/2.0.1-gcc-5.4.0/lib -lmpi
OBJ = acoustics.o acoustics_alg.o export_data.o import_data.o utils.o

all: acoustics
	@mkdir -p export

acoustics: build
	nvcc $(OBJ) $(MPI_LIB) -o acoustics

build:
	nvcc $(MPI_LIB) -c acoustics.cu -o acoustics.o
	nvcc $(MPI_LIB) -c acoustics_alg.cu -o acoustics_alg.o
	nvcc $(MPI_LIB) -c export_data.cu -o export_data.o
	nvcc $(MPI_LIB) -c import_data.cu -o import_data.o
	nvcc $(MPI_LIB) -c utils.cu -o utils.o

clean:
	rm -f *~ *.o acoustics
	rm -fr export
