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


library: includes $(OBJS)
	$(AR) $(ARFLAGS) lib$(SOLVERLIB).a $(OBJS)
	mv lib$(SOLVERLIB).a ../lib


includes:
	ln -f Numerics/*.h  ../include/
	ln -f Utilities/*.h ../include/
	ln -f System/*.h    ../include/
	ln -f System/*/*.h  ../include/

clean:
	rm -rf */*.o */*~ *~ *.mod ../lib/*.a System/*/*.o System/*/*~
