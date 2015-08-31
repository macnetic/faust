CXX=g++

MATRIX_DIR  = $(FAUST_ROOT_DIR)/trunk/devcpp/matrix
PALM_DIR    = $(FAUST_ROOT_DIR)/trunk/devcpp/palm4MSA
INIT_MATLAB = $(FAUST_ROOT_DIR)/trunk/devcpp/init_from_matlab
FAUST_CORE = $(FAUST_ROOT_DIR)/trunk/devcpp/faust_core
MATLAB_INTERFACE_DIR = $(FAUST_ROOT_DIR)/trunk/devcpp/interface_matlab



SRC_CPP = $(wildcard $(PALM_DIR)/*.cpp)  
SRC_CPP += $(wildcard $(MATRIX_DIR)/*.cpp)
SRC_CPP += $(wildcard $(FAUST_CORE)/*.cpp)
SRC_CPP += $(wildcard $(MATLAB_INTERFACE_DIR)/tools/*.cpp)

OBJS_MEX_BASE := $(filter-out $(MATRIX_DIR)/faust_timer.o.mex,$(patsubst %.cpp,%.o.mex,$(SRC_CPP)))
OBJS_MEX_BASE := $(filter-out $(PALM_DIR)/sparse_hierarchical_fact.o.mex,$(OBJS_MEX_BASE))
OBJS_MEX_BASE := $(filter-out $(PALM_DIR)/sparsePalm4MSA.o.mex,$(OBJS_MEX_BASE))    
OBJS_MEX := $(OBJS_MEX_BASE) mexHierarchical_fact.o.mex
OBJS_MEX1 := $(OBJS_MEX_BASE) mexLoadFaust.o.mex
OBJS_MEX2 := $(OBJS_MEX_BASE) mexLoadFaustBis.o.mex
OBJS_MEX3 := $(OBJS_MEX_BASE) mexPalm4MSA.o.mex
OBJS_MEX4 := $(OBJS_MEX_BASE) faust_mex.o.mex


CXXFLAGS = -O2  -I$(MATRIX_DIR) -I$(EIGEN_ROOT_DIR) -I$(PALM_DIR)  -DPROX=2   -I$(FAUST_CORE)  -D_GNU_SOURCE -DMATLAB_MEX_FILE -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O3 -DNDEBUG  -I$(MATLAB_ROOT_DIR)/extern/include -I$(MATLAB_INTERFACE_DIR)/tools #-DMX_COMPAT_32



ifneq ($(COMPILATION_THOMAS),)
  # CXXFLAGS += -DFAUST_SINGLE _
   CXXFLAGS += -D__GEMM_WITH_OPENBLAS__ -I$(OPENBLASDIR) 
   CXXFLAGS += -std=c++11
  # CXXFLAGS += -D__PAS_FIXE__ 
  # CXXFLAGS += -D__COMPILE_TIMERS__
else
   CXXFLAGS += -DFAUST_SINGLE  
   CXXFLAGS += -I$(OPENBLASDIR)/include 
   CXXFLAGS += -D__GEMM_WITH_OPENBLAS__ -I$(OPENBLASDIR) 
   CXXFLAGS +=-std=c++0x  
   #CXXFLAGS +=-D__COMPILE_TIMERS__ 
   #CXXFLAGS +=-D__GEMM_WITH_OPENBLAS__
   #CXXFLAGS +=-D__NICO__
   #CXXFLAGS += -D__PAS_FIXE__
endif

 



LDFLAGS_MEX =  -pthread -Wl,--no-undefined -shared -O2 -Wl,--version-script,${MATLAB_ROOT_DIR}/extern/lib/glnxa64/mexFunction.map -L${MATLAB_ROOT_DIR}/bin/glnxa64 -lmx -lmex -lmat -lstdc++ -lm  -lopenblas -L${OPENBLASDIR}/lib



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


EXECUTABLE= mex mex1 mex3 mex4

all: $(EXECUTABLE) 

mex: mexHierarchical_fact.mexa64 

mex1: mexLoadFaust.mexa64



mex3: mexPalm4MSA.mexa64

mex4 : faust_mex.mexa64

 faust_mex.mexa64:$(OBJS_MEX4)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS_MEX) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN) mexLoadFaust.mexa64 executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

 mexPalm4MSA.mexa64:$(OBJS_MEX3)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS_MEX) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN) mexLoadFaust.mexa64 executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'







mexLoadFaust.mexa64:$(OBJS_MEX1)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS_MEX) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN) mexLoadFaust.mexa64 executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'


mexHierarchical_fact.mexa64:$(OBJS_MEX)
	@printf "MEXHIERARCHICAL_FACT"
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS_MEX) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)MexHierarchical_fact executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'


%.o.mex:%.cpp
	@printf '$(FONT)$(CYAN)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n' 'Compiling' '     ' "$@"
	@$(CXX) -c $(CXXFLAGS) $^ -o $@


clean:
	@printf '$(FONT)$(CYAN)Removing all object files and executable$(RESET_FONT)\n'
	rm -f $(EXECUTABLE)   $(OBJS_MEX) $(OBJS_MEX1)  $(OBJS_MEX2) $(OBJS_MEX3) $(OBJS_MEX4) $(MATLAB_INTERFACE_DIR)/*.o.mex
