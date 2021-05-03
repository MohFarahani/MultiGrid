
# Configuration of the executable
TARGET =mgsolve 

INSTALL_PATH = $(PWD)

# Compiler configuration
CXX      = g++
CXXFLAGS = -O3  -Wall -Winline -Wshadow  -ansi 
COMP = -c 
INC = -I$(INSTALL_PATH)
RM = rm -f


OBJ = main.o \
      Grid.o \
      solver.o\
     
    
default: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ) $(INC)

main.o: main.c solver.c
	$(CXX) $(COMP) $(CXXFLAGS) main.c -o $@ $(INC)

Grid.o: Grid.c Grid.h
	$(CXX) $(COMP) $(CXXFLAGS) Grid.c -o $@ $(INC)

solver.o: solver.c solver.h
	$(CXX) $(COMP) $(CXXFLAGS) solver.c -o $@ $(INC)

clean:
	$(RM) $(OBJ) *~

