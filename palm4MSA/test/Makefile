CXX=g++

MATRIX_DIR  = $(FAUST_ROOT_DIR)/trunk/devcpp/matrix
PALM_DIR    = $(FAUST_ROOT_DIR)/trunk/devcpp/palm4MSA
INIT_MATLAB = $(FAUST_ROOT_DIR)/trunk/devcpp/init_from_matlab
FAUST_CORE = $(FAUST_ROOT_DIR)/trunk/devcpp/faust_core



SRC_CPP = $(wildcard $(PALM_DIR)/*.cpp)  #$(filter-out $(PALM_DIR)/hierarchical_fact.cpp,$(wildcard $(PALM_DIR)/*.cpp))    #$(wildcard $(PALM_DIR)/*.cpp)
SRC_CPP += $(wildcard $(INIT_MATLAB)/*.cpp)
SRC_CPP += $(wildcard $(MATRIX_DIR)/*.cpp)
ifneq ($(COMPILATION_THOMAS),)
SRC_CPP += $(wildcard $(FAUST_CORE)/*.cpp)
endif
OBJS := $(patsubst %.cpp,%.o,$(SRC_CPP)) 

CXXFLAGS = -O2  -I$(MATRIX_DIR) -I$(EIGEN_ROOT_DIR) -I$(PALM_DIR) -I$(INIT_MATLAB)  -I$(MATIO_ROOT_DIR)/include  -DPROX=2   -I$(FAUST_CORE) #-Winline 

ifneq ($(COMPILATION_THOMAS),)
  # CXXFLAGS += -DFAUST_SINGLE _
   CXXFLAGS += -D__GEMM_WITH_OPENBLAS__ -I$(OPENBLASDIR) 
   CXXFLAGS += -std=c++11
  # CXXFLAGS += -D__PAS_FIXE__ 
  # CXXFLAGS += -D__COMPILE_TIMERS__
else
   CXXFLAGS += -DFAUST_SINGLE  
   CXXFLAGS += -I$(OPENBLASDIR)/include 
   CXXFLAGS +=-std=c++0x  
   #CXXFLAGS +=-D__COMPILE_TIMERS__ 
   CXXFLAGS +=-D__GEMM_WITH_OPENBLAS__
   CXXFLAGS +=-D__NICO__
   CXXFLAGS += -D__PAS_FIXE__
endif

LDFLAGS = -L$(MATIO_ROOT_DIR)/lib -lmatio 
ifneq ($(COMPILATION_THOMAS),)
   LDFLAGS += -lhdf5  -L$(OPENBLASDIR) -lopenblas
else
   LDFLAGS += -lrt -L$(OPENBLASDIR)/lib   -lopenblas
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


all:hier_fact hier_fact2 palm meg meg_new meg_old faust

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



%.o:%.cpp
	@printf '$(FONT)$(CYAN)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n' 'Compiling' '     ' "$@"
	@$(CXX) -c $(CXXFLAGS) $^ -o $@

clean:
	@printf '$(FONT)$(CYAN)Removing all object files and executable$(RESET_FONT)\n'
	@rm -f hier_fact hier_fact2 palm meg meg_new meg_old faust  $(MATRIX_DIR)/*.o $(PALM_DIR)/*.o $(INIT_MATLAB)/*.o $(FAUST_CORE)/*.o ./*.o
