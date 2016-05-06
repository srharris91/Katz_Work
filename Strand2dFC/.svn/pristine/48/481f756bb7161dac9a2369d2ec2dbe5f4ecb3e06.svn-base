###############################################################################
#                                                                             #
#                    Top level Makefile for arccfd codes                      #
#                                                                             #
###############################################################################


library:
	@echo ""
	@echo "making library . . ."
	@echo ""
	cd src; make library


documentation:
	@echo ""
	@echo "making documentation . . ."
	@echo ""
	cd doc; doxygen Doxyfile


clean: cleanCode cleanDoc
	rm -f *~ */*~


cleanCode:
	@echo ""
	@echo "cleaning code . . ."
	@echo ""
	cd src; make clean
	cd lib; rm -rf *.a *~
	cd include; rm -rf *.h *~


cleanDoc:
	@echo ""
	@echo "cleaning documentation . . ."
	@echo ""
	cd doc; rm -rf html latex *~
