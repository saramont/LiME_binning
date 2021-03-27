VLIB= -g -O0

MY_CXX_FLAGS= -std=c++11 -Wall -DNDEBUG -fomit-frame-pointer -Wno-comment

OMP_LIB= -fopenmp 

MY_CXX_OPT_FLAGS= -O3 -m64 $(OMP_LIB) 

MY_CXX=g++


OMP = 0
DEBUG = 0
POS = 0


LFLAGS = -lm -ldl

DEFINES = -DDEBUG=$(DEBUG) -DOMP=$(OMP) -DPOS=$(POS)

CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(LFLAGS) $(DEFINES) -I$(INC_DIR) -L$(LIB_DIR)

##

all: compile

clean:
	\rm -f *.o ClusterLCP BinningDA modifyDA

##

compile: ClusterLCP BinningDA modifyDA

ClusterLCP: src/ClusterLCP.cpp ${LIBOBJ} 
	$(MY_CXX) src/ClusterLCP.cpp $(CCLIB) -o ClusterLCP ${LIBOBJ} $(CXX_FLAGS) 

BinningDA: src/BinningDA.cpp ${LIBOBJ} 
	$(MY_CXX) src/BinningDA.cpp $(CCLIB) -o BinningDA ${LIBOBJ} $(CXX_FLAGS) 

modifyDA: src/modifyDA.cpp ${LIBOBJ} 
	$(MY_CXX) src/modifyDA.cpp $(CCLIB) -o modifyDA ${LIBOBJ} $(CXX_FLAGS) 
