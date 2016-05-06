###############################################################################
#                                                                             #
#                    Src level Makefile for arccfd codes                      #
#                                                                             #
###############################################################################


include ../Make.rules
include Numerics/Make.Numerics.objs
include System/Make.System.objs
include Utilities/Make.Utilities.objs


OBJS = $(NUMERICSOBJS) $(SYSTEMOBJS) $(UTILITIESOBJS)


driver: library ../run/main.o
	$(LINKER) $(CXXFLAGS) $(CXXFLAGS_EXTRA) $(INCLUDE) \
	../run/main.o -o $(EXEC) \
	-L$(MESHER)/lib  -l$(MESHERLIB)  \
	-L../lib         -l$(SOLVERLIB)  \
	-L$(PHYSICS)/lib -l$(PHYSICSLIB) \
	$(LDFLAGS)
	mv $(EXEC) ../run


includes:
	ln -f Numerics/*.h  Include/
	ln -f Utilities/*.h Include/
	ln -f System/*.h    Include/
	ln -f System/*/*.h  Include/


library: includes $(OBJS)
	$(AR) $(ARFLAGS) lib$(SOLVERLIB).a $(OBJS)
	mv lib$(SOLVERLIB).a ../lib


clean:
	rm -rf ../run/$(EXEC) ../run/*~ *.dSYM \
	*/*.o */*~ *~ *.mod ../lib/*.a System/*/*.o System/*/*~