CXX=g++

MATRIX_DIR  = $(FAUST_ROOT_DIR)/trunk/devcpp/matrix
PALM_DIR    = $(FAUST_ROOT_DIR)/trunk/devcpp/palm4MSA
INIT_MATLAB = $(FAUST_ROOT_DIR)/trunk/devcpp/init_from_matlab
FAUST_CORE = $(FAUST_ROOT_DIR)/trunk/devcpp/faust_core
MATLAB_INTERFACE_DIR = $(FAUST_ROOT_DIR)/trunk/devcpp/interface_matlab
MATLAB_INTERFACE_TOOL_DIR = $(MATLAB_INTERFACE_DIR)/tools



SRC_CPP = $(wildcard $(PALM_DIR)/*.cpp)  #$(filter-out $(PALM_DIR)/hierarchical_fact.cpp,$(wildcard $(PALM_DIR)/*.cpp))    #$(wildcard $(PALM_DIR)/*.cpp)
SRC_CPP += $(wildcard $(INIT_MATLAB)/*.cpp)
SRC_CPP += $(wildcard $(MATRIX_DIR)/*.cpp)
#ifneq ($(COMPILATION_THOMAS),)
SRC_CPP += $(wildcard $(FAUST_CORE)/*.cpp)
#endif
OBJS := $(patsubst %.cpp,%.o,$(SRC_CPP)) 
OBJS_MEX_COMMON := $(filter-out $(MATRIX_DIR)/faust_timer.cpp,$(patsubst %.cpp,%.o.mex,$(SRC_CPP))) $(MATLAB_INTERFACE_TOOL_DIR)/mexFaustMat.o.mex 


CXXFLAGS = -O2  -I$(MATRIX_DIR) -I$(EIGEN_ROOT_DIR) -I$(PALM_DIR) -I$(INIT_MATLAB)  -I$(MATIO_ROOT_DIR)/include  -DPROX=2   -I$(FAUST_CORE) #-Winline 

CXXFLAGS_MEX = $(CXXFLAGS) -DMX_COMPAT_32 -D_GNU_SOURCE -DMATLAB_MEX_FILE -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O3 -DNDEBUG  -I$(MATLAB_ROOT_DIR)/extern/include -I$(MATLAB_INTERFACE_TOOL_DIR)

ifneq ($(COMPILATION_THOMAS),)
   CXXFLAGS += -DFAUST_SINGLE 
   CXXFLAGS += -D__GEMM_WITH_OPENBLAS__ -I$(OPENBLASDIR) 
   CXXFLAGS += -std=c++11
  # CXXFLAGS += -D__PAS_FIXE__ 
  # CXXFLAGS += -D__COMPILE_TIMERS__
else
   CXXFLAGS += -DFAUST_SINGLE  
   #CXXFLAGS += -I$(OPENBLASDIR)/include 
   CXXFLAGS +=-std=c++0x  
   CXXFLAGS +=-D__COMPILE_TIMERS__ 
   #CXXFLAGS +=-D__GEMM_WITH_OPENBLAS__
   CXXFLAGS +=-D__NICO__
   #CXXFLAGS += -D__PAS_FIXE__
endif

LDFLAGS = -L$(MATIO_ROOT_DIR)/lib -lmatio 
LDFLAGS_MEX = $(LDFLAGS) -pthread -Wl,--no-undefined -shared -O2 -Wl,--version-script,${MATLAB_ROOT_DIR}/extern/lib/glnxa64/mexFunction.map -L${MATLAB_ROOT_DIR}/bin/glnxa64 -lmx -lmex -lmat -lstdc++ -lm 

ifneq ($(COMPILATION_THOMAS),)
   LDFLAGS += -lhdf5  -L$(OPENBLASDIR) -lopenblas
else
   LDFLAGS += -lhdf5 -L$(HDF5_ROOT_DIR)/lib -lrt -L$(OPENBLASDIR)/lib   -lopenblas
endif

FONT =\033[

BOLD      =1;
ITALIC    =3;
UNDERLINE =4;
BLINK     =5;

RED     =31m
GREEN   =32m
YELLOW  =33m #orange in non-bold or yellow in bold
BLUE    =34m
MAGENTA =35m
CYAN    =36m

RESET_FONT =\033[0m


MEX_EXE = $(MATLAB_INTERFACE_DIR)/mexHierarchical_fact.mexa64 $(MATLAB_INTERFACE_DIR)/mexLoadFaust.mexa64 $(MATLAB_INTERFACE_DIR)/faust_mex.mexa64

EXECUTABLE=hier_fact hier_fact2 palm meg meg_new meg_old faust comptime0 comptime comptime2 $(MEX_EXE) 

# instruction pour que les fichiers .o.mex ne soient pas automatiquement detruits lors apres les appels aux cibles all ou mex
.PRECIOUS: %.o.mex

all: $(EXECUTABLE)

mex: $(MEX_EXE)

hier_fact: hierarchical_fact_test.o $(OBJS)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

hier_fact2: hierarchical_fact_test2.o $(OBJS)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'



palm:$(OBJS) palm4MSA_test.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

meg:$(OBJS) MEG_fact.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'
	
meg_new:$(OBJS) MEG_fact.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'	
	
meg_old:$(OBJS) MEG_fact.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'		

faust:$(OBJS) faust_test.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

comptime0:$(OBJS) comp_time0.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'
	
comptime:$(OBJS) comp_time.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'
	
comptime2:$(OBJS) comp_time2.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

%.mexa64: %.o.mex $(OBJS_MEX_COMMON) 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS_MEX) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Mex executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'




%.o:%.cpp
	@printf '$(FONT)$(CYAN)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n' 'Compiling' '     ' "$@"
	@$(CXX) -c $(CXXFLAGS) $^ -o $@

%.o.mex:%.cpp
	@printf '$(FONT)$(CYAN)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n' 'Compiling' '     ' "$@"
	@$(CXX) -c $(CXXFLAGS_MEX) $^ -o $@


clean:
	@printf '$(FONT)$(CYAN)Removing all object files and executable$(RESET_FONT)\n'
	@rm -f $(EXECUTABLE)  $(OBJS) $(OBJS_MEX_COMMON) *.o $(MATLAB_INTERFACE_DIR)/*.o.mex
