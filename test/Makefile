CXX=g++

MATRIX_DIR  = /home/tgautrai/faust2/matrix
PALM_DIR    = /home/tgautrai/faust2/palm4MSA
INIT_MATLAB = /home/tgautrai/faust2/init_from_matlab
FAUST_CORE = /home/tgautrai/faust2/faust_core
ERROR_DIR = /home/tgautrai/faust2/exception
TESTS_DIR = /home/tgautrai/faust2/test/src
MATLAB_INTERFACE_DIR = /home/tgautrai/faust2/interface_matlab
MATLAB_INTERFACE_TOOL_DIR = $(MATLAB_INTERFACE_DIR)/tools
MATLAB_INTERFACE_MEX_DIR = $(MATLAB_INTERFACE_DIR)/mex_functions

SRC_CPP = $(wildcard $(PALM_DIR)/*.cpp)  #$(filter-out $(PALM_DIR)/hierarchical_fact.cpp,$(wildcard $(PALM_DIR)/*.cpp))    #$(wildcard $(PALM_DIR)/*.cpp)
SRC_CPP += $(wildcard $(INIT_MATLAB)/*.cpp)
SRC_CPP += $(wildcard $(MATRIX_DIR)/*.cpp)
#ifneq ($(COMPILATION_THOMAS),)
SRC_CPP += $(wildcard $(FAUST_CORE)/*.cpp)
#endif
SRC_CPP += $(wildcard $(ERROR_DIR)/*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRC_CPP)) 
OBJS_MEX_COMMON := $(filter-out $(MATRIX_DIR)/faust_timer.cpp,$(patsubst %.cpp,%.o.mex,$(SRC_CPP))) $(MATLAB_INTERFACE_TOOL_DIR)/mexFaustMat.o.mex 


CXXFLAGS = -O2 -Wall -I$(MATRIX_DIR) -I$(EIGEN_ROOT_DIR) -I$(PALM_DIR) -I$(INIT_MATLAB)  -I$(MATIO_ROOT_DIR)/include  -DPROX=2   -I$(FAUST_CORE) -I$(ERROR_DIR) #-Winline  

CXXFLAGS_MEX = $(CXXFLAGS) -D_GNU_SOURCE -DMATLAB_MEX_FILE -DCOMPILE_MEX -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O3 -DNDEBUG  -I$(MATLAB_ROOT_DIR)/extern/include -I$(MATLAB_INTERFACE_MEX_DIR) -I$(MATLAB_INTERFACE_TOOL_DIR) -I$(ERROR_DIR)

ifneq ($(COMPILATION_THOMAS),)
   CXXFLAGS += -DFAUST_SINGLE 
   CXXFLAGS += -D__GEMM_WITH_OPENBLAS__ -I$(OPENBLASDIR) 
   CXXFLAGS += -std=c++11
  # CXXFLAGS += -D__PAS_FIXE__ 
  # CXXFLAGS += -D__COMPILE_TIMERS__
else
   #CXXFLAGS += -DFAUST_SINGLE  
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
   LDFLAGS += -L$(HDF5_ROOT_DIR)/lib -lrt -L$(OPENBLASDIR)/lib   -lopenblas  -lhdf5
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


MEX_EXE = $(MATLAB_INTERFACE_MEX_DIR)/mexHierarchical_fact.mexa64 $(MATLAB_INTERFACE_MEX_DIR)/mexLoadFaust.mexa64 $(MATLAB_INTERFACE_MEX_DIR)/faust_mex.mexa64

#EXECUTABLE=$(TESTS_DIR)/hier_fact $(TESTS_DIR)/hier_fact2 $(TESTS_DIR)/hier_fact3 $(TESTS_DIR)/palm $(TESTS_DIR)/meg $(TESTS_DIR)/meg_new $(TESTS_DIR)/meg_old $(TESTS_DIR)/faust $(TESTS_DIR)/comptime0 $(TESTS_DIR)/comptime $(TESTS_DIR)/comptime2 $(MEX_EXE) 
EXECUTABLE=$(TESTS_DIR)/comptime0 $(TESTS_DIR)/comptime $(TESTS_DIR)/comptime2 

# instruction pour que les fichiers .o.mex ne soient pas automatiquement detruits apres les appels aux cibles all ou mex
.PRECIOUS: %.o.mex

all: regle_bidon $(EXECUTABLE) 

regle_bidon :
	@printf "$(SRC_CPP)"

mex: $(MEX_EXE)

$(TESTS_DIR)/hier_fact: $(TESTS_DIR)/hierarchical_fact_test.o $(OBJS)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

$(TESTS_DIR)/hier_fact2: $(TESTS_DIR)/hierarchical_fact_test2.o $(OBJS)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'
	
$(TESTS_DIR)/hier_fact3: $(TESTS_DIR)/hierarchical_fact_test3.o $(OBJS)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'



$(TESTS_DIR)/palm:$(OBJS) $(TESTS_DIR)/palm4MSA_test.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

$(TESTS_DIR)/meg:$(OBJS) $(TESTS_DIR)/MEG_fact.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'
	
$(TESTS_DIR)/meg_new:$(OBJS) $(TESTS_DIR)/MEG_fact.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'	
	
$(TESTS_DIR)/meg_old:$(OBJS) $(TESTS_DIR)/MEG_fact.o
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'		

$(TESTS_DIR)/faust:$(OBJS) $(TESTS_DIR)/faust_test.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

$(TESTS_DIR)/comptime0:$(OBJS) $(TESTS_DIR)/comp_time0.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'
	
$(TESTS_DIR)/comptime:$(OBJS) $(TESTS_DIR)/comp_time.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'
	
$(TESTS_DIR)/comptime2:$(OBJS) $(TESTS_DIR)/comp_time2.o 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

%.mexa64: %.o.mex $(OBJS_MEX_COMMON) 
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS_MEX) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Mex executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'




%.o: %.cpp
	@printf 'CXXFLAGS $(CXXFLAGS)'
	@printf '$(FONT)$(CYAN)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n' 'Compiling' '     ' "$@"
	@$(CXX) -c $(CXXFLAGS) $^ -o $@

%.o.mex: %.cpp
	@printf '$(FONT)$(CYAN)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n' 'Compiling' '     ' "$@"
	@$(CXX) -c $(CXXFLAGS_MEX) $^ -o $@


clean:
	@printf '$(FONT)$(CYAN)Removing all object files and executable$(RESET_FONT)\n'
	@rm -f $(EXECUTABLE)  $(OBJS) $(OBJS_MEX_COMMON) *.o $(MATLAB_INTERFACE_MEX_DIR)/*.o.mex $(MATLAB_INTERFACE_TOOL_DIR)/*.o.mex $(ERROR_DIR)/*.o.mex
