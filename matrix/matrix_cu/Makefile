CC = gcc
NVCC = nvcc

CFLAGS = -O3 -fopenmp -I$(EIGENDIR)  -m64 -DHAVE_CUBLAS -fPIC -I$(CUDADIR)/include -fopenmp

CUFLAGS = -O3 -m64 -Xcompiler -fPIC -use_fast_math -arch=compute_20 -code=sm_21 -I$(CUDADIR)/include -DHAVE_CUBLAS

LIB_FAUST_CLASSES = libfaust_classes.a



SRC_CPP := $(wildcard ./*.cpp)
SRC_CU  := $(wildcard ./*.cu)
SRC_HPP := $(wildcard ./*.hpp)
SRC_HU  := $(wildcard ./*.hu)

ifneq ($(COMPILE_GPU),)
  CUFLAGS += -D__COMPILE_GPU__
  CFLAGS  += -D__COMPILE_GPU__
endif

ifneq ($(COMPILE_SPMAT),)
  CUFLAGS += -D__COMPILE_SPMAT__
  CFLAGS  += -D__COMPILE_SPMAT__
else
  SRC_CPP := $(filter-out ./faust_spmat.cpp ./faust_cu_spmat.cpp,$(SRC_CPP))
endif

OBJ_CU  := $(patsubst ./%.cu, ./%.cuo, $(SRC_CU))
OBJ_CPP := $(patsubst ./%.cpp, ./%.o, $(SRC_CPP))


LDFLAGS =  -lpthread -lm -lcublas -lcudart -lcusparse -lstdc++ -lgfortran -lz -lmatio


all: $(LIB_FAUST_CLASSES)

./%.cuo:  ./%.cu $(SRC_HU)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(CUFLAGS) $(LDFLAGS) -c -o $@ $<
./%.o: ./%.cpp $(SRC_HPP)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) $(CFLAGS) $(LDFLAGS) -c -o $@ $< 

$(LIB_FAUST_CLASSES): $(OBJ_CPP) $(OBJ_CU) $(SRC_HPP)
	@ar -rv $@ $(OBJ_CPP) $(OBJ_CU)


clean:
	@rm -f $(OBJ_CU) $(OBJ_CPP) $(LIB_FAUST_CLASSES)


