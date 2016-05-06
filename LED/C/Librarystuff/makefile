#~ executable:
	#~ mkdir ./EXECUTABLE
	#~ g++ librarycode.cpp -o librarycode.exec
	#~ mv librarycode.exec ./EXECUTABLE
	#~ cd ./EXECUTABLE && ./librarycode.exec
#~ clean:
	#~ $(RM) -r EXECUTABLE/


## simple make file  
CC=g++ #compiler option
CFLAGS=-c -Wall#flags for compiler
LDFLAGS= 
SOURCES=librarycode.cpp librarycode2.cpp
OBJECTS=$(SOURCES:.cpp=.o) #make object files of sources names ($ is dereference operator)
EXECUTABLE=librarycode.exec #single executable file name

all: $(SOURCES) $(EXECUTABLE)
	./$(EXECUTABLE) #execute program
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE) #link to make executable
#~ -o $@ says to put the output of the compilation in the file named on the left side of the :
#~ $^ references things on the right of :
library: $(OBJECTS)
	ar cr libname.a $(OBJECTS) && $(CC) $(CFLAGS) $(OBJECTS) \
	 -L. -llibname.a -o $(EXECUTABLE)
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ #compile to make object files
# $< is first item of the list

clean:
	rm -rf *.exec *.o *.a #remove .exec .o and .a files 
