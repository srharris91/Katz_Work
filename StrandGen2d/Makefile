###############################################################################
#                                                                             #
#                    Src level Makefile for codes                             #
#                                                                             #
###############################################################################


include Make.rules

STRANDGEN_OBJS = mesh/*.o		\
	         solver/*.o

driver: library driver/main.o
	$(LINKER) $(CXXFLAGS) $(INCLUDE) driver/main.o \
	-o $(EXEC) -L./lib -l$(STRANDGENLIB) $(LDFLAGS)
	mv $(EXEC) ./driver/
	rm -rf *.dSYM


library:
	( cd mesh    && $(MAKE) all ) || exit 1
	( cd solver  && $(MAKE) all ) || exit 1
#	( cd utils   && $(MAKE) all ) || exit 1
	$(AR) r lib$(STRANDGENLIB).a $(STRANDGEN_OBJS)
	mv lib$(STRANDGENLIB).a lib/


clean:
	( cd mesh    && $(MAKE) clean ) || exit 1
	( cd solver  && $(MAKE) clean ) || exit 1
#	( cd utils   && $(MAKE) clean ) || exit 1
	$(RM) lib/lib$(STRANDGENLIB).a
	$(RM) *.o */*.o *~ */*~ *.mod */*.mod

