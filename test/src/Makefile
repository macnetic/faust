CC = gcc
NVCC = nvcc


MATRIX_DIR = ../../matrix/matrix_cu
PALM_DIR = ../../palm4MSA/palm4MSA_cu
FAUSTCORE_DIR = ../../faust_core/faust_core_cu
MATIO_DIR = ../../init_from_matlab
EXCEPTION_DIR = ../../exception
TEST_DIR = .

LIB_MATRIX = $(MATRIX_DIR)/libfaustmat.a
LIB_PALM = $(PALM_DIR)/libfaustpalm.a
LIB_FAUSTCORE = $(FAUSTCORE_DIR)/libfaustcore.a
LIB_MATIO = $(MATIO_DIR)/libfaustmatio.a
LIB_EXCEPTION = $(EXCEPTION_DIR)/libfaustexeption.a

LIBS = $(LIB_MATRIX) $(LIB_PALM) $(LIB_FAUSTCORE) $(LIB_MATIO) $(LIB_EXCEPTION)

CFLAGS = -O3 -fopenmp -m64 -DHAVE_CUBLAS -fPIC
CUFLAGS = -O3 -m64 -Xcompiler -fPIC -use_fast_math -arch=compute_20 -code=sm_21

ifneq ($(COMPILE_GPU),)
  CUFLAGS += -D__COMPILE_GPU__
  CFLAGS  += -D__COMPILE_GPU__
endif
ifneq ($(COMPILE_SPMAT),)
  CUFLAGS += -D__COMPILE_SPMAT__
  CFLAGS  += -D__COMPILE_SPMAT__
endif


CUFLAGS_MATRIX   := $(CUFLAGS) 

CFLAGS_MATRIX    := $(CFLAGS) -I$(CUDADIR)/include -I$(EXCEPTION_DIR) 
CFLAGS_PALM      := $(CFLAGS) -I$(MATRIX_DIR) -I$(CUDADIR)/include -I$(EIGENDIR) -I$(EXCEPTION_DIR)
CFLAGS_FAUSTCORE := $(CFLAGS) 
CFLAGS_MATIO     := $(CFLAGS) -I$(MATIODIR)/include -I$(PALM_DIR) -I$(MATRIX_DIR) -I$(CUDADIR)/include -I$(EIGENDIR) -I$(EXCEPTION_DIR)
CFLAGS_EXCEPTION := $(CFLAGS) 
CFLAGS_TEST      := $(CFLAGS) -I$(MATRIX_DIR) -I$(CUDADIR)/include -I$(EIGENDIR) -I$(PALM_DIR) -I$(FAUSTCORE_DIR) -I$(MATIO_DIR) -I$(MATIODIR)/include 

SRC_CPP_MATRIX := $(wildcard $(MATRIX_DIR)/*.cpp)
SRC_CU_MATRIX  := $(wildcard $(MATRIX_DIR)/*.cu)
SRC_HPP_MATRIX := $(wildcard $(MATRIX_DIR)/*.hpp)
SRC_HU_MATRIX  := $(wildcard $(MATRIX_DIR)/*.hu)

SRC_CPP_PALM := $(wildcard $(PALM_DIR)/*.cpp)
SRC_HPP_PALM := $(wildcard $(PALM_DIR)/*.hpp)

SRC_CPP_FAUSTCORE := $(wildcard $(LIB_FAUSTCORE)/*.cpp)
SRC_HPP_FAUSTCORE := $(wildcard $(LIB_FAUSTCORE)/*.hpp)

SRC_CPP_MATIO := $(wildcard $(MATIO_DIR)/*.cpp)
SRC_HPP_MATIO := $(wildcard $(MATIO_DIR)/*.hpp)

SRC_CPP_EXCEPTION := $(wildcard $(EXCEPTION_DIR)/*.cpp)
SRC_HPP_EXCEPTION := $(wildcard $(EXCEPTION_DIR)/*.hpp)


ifneq ($(COMPILE_SPMAT),)
  SRC_CPP_MATRIX := $(filter-out $(MATRIX_DIR)/faust_spmat.cpp $(MATRIX_DIR)/faust_cu_spmat.cpp,$(SRC_CPP_MATRIX))
endif
SRC_CPP_MATRIX := $(filter-out $(MATRIX_DIR)/cblas_algebra.cpp,$(SRC_CPP_MATRIX))


OBJ_CU_MATRIX     := $(patsubst %.cu,  %.cuo, $(SRC_CU_MATRIX))

OBJ_CPP_MATRIX    := $(patsubst %.cpp, %.o, $(SRC_CPP_MATRIX))
OBJ_CPP_PALM      := $(patsubst %.cpp, %.o, $(SRC_CPP_PALM))
OBJ_CPP_FAUSTCORE := $(patsubst %.cpp, %.o, $(SRC_CPP_FAUSTCORE))
OBJ_CPP_MATIO     := $(patsubst %.cpp, %.o, $(SRC_CPP_MATIO))
OBJ_CPP_EXCEPTION := $(patsubst %.cpp, %.o, $(SRC_CPP_EXCEPTION))


TEST_HIERFACT_SRC = hierarchical_fact_test_cu.cpp
TEST_HIERFACT_OBJ = $(patsubst %.cpp, %.o, $(TEST_HIERFACT_SRC))
TEST_HIERFACT_EXE = $(patsubst %.o, %.out, $(TEST_HIERFACT_OBJ))

LDFLAGS = $(LIBS) -L$(CUDADIR)/lib64 -L$(MATIODIR)/lib -lpthread -lm -lcublas -lcudart -lcusparse -lstdc++ -lgfortran -lz -lmatio -lhdf5


all: $(TEST_HIERFACT_EXE)
$(TEST_HIERFACT_EXE) : $(TEST_HIERFACT_OBJ)
$(TEST_HIERFACT_OBJ) : $(TEST_HIERFACT_SRC)



%.out: %.o $(LIB_MATRIX) $(LIB_PALM) $(LIB_FAUSTCORE) $(LIB_MATIO) $(LIB_EXCEPTION)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) $< $(CFLAGS) $(LDFLAGS) -o $@ 


$(MATRIX_DIR)/%.cuo:  $(MATRIX_DIR)/%.cu $(SRC_HU_MATRIX)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) -c $(CUFLAGS_MATRIX) $< -o $@
$(MATRIX_DIR)/%.o: $(MATRIX_DIR)/%.cpp $(SRC_HPP_MATRIX)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) -c $(CFLAGS_MATRIX) $< -o $@
$(LIB_MATRIX): $(OBJ_CPP_MATRIX) $(OBJ_CU_MATRIX) $(SRC_HPP_MATRIX)
	@ar -rv $@ $(OBJ_CPP_MATRIX) $(OBJ_CU_MATRIX)

$(PALM_DIR)/%.o: $(PALM_DIR)/%.cpp $(SRC_HPP_PALM)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) -c $(CFLAGS_PALM) $< -o $@
$(LIB_PALM): $(OBJ_CPP_PALM) $(SRC_HPP_PALM)
	@ar -rv $@ $(OBJ_CPP_PALM)

$(FAUSTCORE_DIR)/%.o: $(FAUSTCORE_DIR)/%.cpp $(SRC_HPP_FAUSTCORE)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) -c $(CFLAGS_FAUSTCORE) $< -o $@
$(LIB_FAUSTCORE): $(OBJ_CPP_FAUSTCORE) $(SRC_HPP_FAUSTCORE)
	@ar -rv $@ $(OBJ_CPP_FAUSTCORE)

$(MATIO_DIR)/%.o: $(MATIO_DIR)/%.cpp $(SRC_HPP_MATIO)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) -c $(CFLAGS_MATIO) $< -o $@
$(LIB_MATIO): $(OBJ_CPP_MATIO) $(SRC_HPP_MATIO)
	@ar -rv $@ $(OBJ_CPP_MATIO)

$(EXCEPTION_DIR)/%.o: $(EXCEPTION_DIR)/%.cpp $(SRC_HPP_EXCEPTION)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) -c $(CFLAGS_EXCEPTION) $< -o $@
$(LIB_EXCEPTION): $(OBJ_CPP_EXCEPTION) $(SRC_HPP_EXCEPTION)
	@ar -rv $@ $(OBJ_CPP_EXCEPTION)

$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) -c $(CFLAGS_TEST) $< -o $@



cleanall:
	@rm -rf  $(LIB_MATRIX) $(LIB_PALM) $(LIB_FAUSTCORE) $(LIB_MATIO) $(LIB_EXCEPTION) $(MATRIX_DIR)/*.o $(MATRIX_DIR)/*.cuo $(PALM_DIR)/*.o $(FAUSTCORE_DIR)/*.o $(MATIO_DIR)/*.o $(EXCEPTION_DIR)/*.o

clean:
	@rm -rf $(TEST_HIERFACT_EXE) $(TEST_HIERFACT_OBJ) 
