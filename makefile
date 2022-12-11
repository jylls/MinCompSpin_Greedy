########################################################################################################################
####################################      CONSTANT TO SPECIFY     ######################################################
########################################################################################################################
#### If the variables are left empty:           ########################################################################
####     then the program will automatically use the filename and number of variables specifed in data.h    ############
########################################################################################################################

datafilename := INPUT/SCOTUS_n9_N895_Data.dat  # datafile name
n := 9		# number of binary variables in the datafile 


########################################################################################################################
#### TO COMPILE:  	Enter "make" in your terminal
#### TO RUN: 		Enter "make run" in your terminal
########################################################################################################################


########################################################################################################################
##########################################      DO NOT MODIFY     ######################################################
########################################################################################################################

CC = g++ 	# Flag for implicit rules: used for linker
CXX = g++ 	# Flag for implicit rules: compilation of c++ files
CXXFLAGS = -std=c++11 -O2  #Extra flags to give to the C++ compiler

# !!!! Implicit compilation of the all the objects !!!!
# Compilation -- Implicite rule:   
#%.o : %.c
#		$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

directory = MCM
objects = tools.o LogE.o LogL.o Complexity.o MCM_info.o Basis_Choice.o P_s.o Operations_OnData.o info_quant.o BestMCM_GreedySearch.o

# String substitution: add the directory name:
# As an example, hello.o turns into ./MCM/hello.o
OBJS := $(objects:%=$(directory)/%)

a.out: $(OBJS) main.o main_routines.o library.hpp.gch
	g++ -std=c++11 main.o main_routines.o $(OBJS)

main.o: main.cpp library.hpp.gch data.h
	g++ -std=c++11 -c main.cpp -include library.hpp # Compile main.cpp

main_routines.o: main_routines.cpp library.hpp.gch
	g++ -std=c++11 -c main_routines.cpp -include library.hpp # Compile main_routines.cpp

library.hpp.gch: library.hpp
	g++ -std=c++11 -c library.hpp

run: a.out
	./a.out $(datafilename) $n

clean:
	rm $(OBJS) library.hpp.gch main_routines.o main.o a.out
