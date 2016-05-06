include ../Make.rules

DIR   = mesh

CXX_OBJS  = StrandBlock.o		\
	    StrandGlobalMesh.o		\
	    StrandMultiBlockMesh.o

F90_OBJS  = avdefs.o			\
	    readinput.o			\
	    readmeshheader.o		\
	    readmeshdata.o		\
	    strand1dDist.o		\
	    sfuns.o			\
	    partitionstage1.o		\
	    partitionstage2.o		\
	    findcspdim.o		\
	    fillcsp.o			\
	    findpspdim.o		\
	    fillpsp.o			\
	    smoothingiter.o		\
	    strandclipsimple.o		\
	    plotheader.o		\
	    plotblock.o			\
	    findsharpcorners.o		\
	    fillsharpcorners.o


all:	$(CXX_OBJS) $(F90_OBJS) 

clean:
	$(RM) *.o *~

