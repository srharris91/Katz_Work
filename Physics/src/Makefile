###############################################################################
#                                                                             #
#                    Src level Makefile for arccfd codes                      #
#                                                                             #
###############################################################################


include ../Make.rules
include State/Make.State.objs
include Transport/Make.Transport.objs
include Solution/Make.Solution.objs
include Fortran/Make.Fortran.objs


OBJS = $(STATEOBJS) $(TRANSPORTOBJS) $(SOLUTIONOBJS) $(FORTRANOBJS)


includes:
	rm -rf ../include
	mkdir ../include
	ln -f *.h           ../include/
	ln -f State/*.h     ../include/
	ln -f Transport/*.h ../include/
	ln -f Solution/*.h  ../include/


library: includes $(OBJS)
	$(AR) $(ARFLAGS) lib$(PHYSICSLIB).a $(OBJS)
	mv lib$(PHYSICSLIB).a ../lib


clean:
	rm -rf */*.o */*~ *~ ../lib/*.a
