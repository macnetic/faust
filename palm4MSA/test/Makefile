CXX=g++

MATRIX_DIR = ../../matrix
PALM_DIR = ..
INIT_MATLAB = ../../init_from_matlab

SRC_CPP = hierarchical_fact_test.cpp    #palm4MSA_test.cpp 
SRC_CPP += $(wildcard $(PALM_DIR)/*.cpp)  #$(filter-out $(PALM_DIR)/hierarchical_fact.cpp,$(wildcard $(PALM_DIR)/*.cpp))    #$(wildcard $(PALM_DIR)/*.cpp)
SRC_CPP += $(wildcard $(INIT_MATLAB)/*.cpp)
SRC_CPP += $(wildcard $(MATRIX_DIR)/*.cpp)

OBJS := $(patsubst %.cpp,%.o,$(SRC_CPP)) 

CXXFLAGS = -O3  -I$(MATRIX_DIR) -I$(EIGENDIR) -I$(PALM_DIR) -I$(INIT_MATLAB) -I$(MATIODIR)/include #-Winline

LDFLAGS = -L$(MATIODIR)/lib -lmatio -lhdf5

FONT =\033[

BOLD      =1;
ITALIC    =3;
UNDERLINE =4;
BLINK     = 5;

RED     =31m
GREEN   =32m
YELLOW  =33m#orange in non-bold or yellow in bold
BLUE    =34m
MAGENTA =35m
CYAN    =36m

RESET_FONT =\033[0m




all:test

test:$(OBJS)
	@printf "$(FONT)$(YELLOW)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n" 'Linking' '       ' "$@"
	@$(CXX) $(LDFLAGS) $^ -o $@
	@printf '\n$(FONT)$(BOLD)$(GREEN)Executable has been correctly generated$(FONT)$(BOLD)$(BLINK)$(YELLOW) :)$(RESET_FONT)\n\n\n'

%.o:%.cpp
	@printf '$(FONT)$(CYAN)%s$(RESET_FONT)%s$(FONT)$(ITALIC)$(RED)%s$(RESET_FONT)\n' 'Compiling' '     ' "$@"
	@$(CXX) -c $(CXXFLAGS) $^ -o $@

clean:
	@printf '$(FONT)$(CYAN)Removing all object files and executable$(RESET_FONT)\n'
	@rm -f test $(MATRIX_DIR)/*.o $(PALM_DIR)/*.o $(INIT_MATLAB)/*.o ./*.o
