###############################################################################
#                                                                             #
#                    Top level Makefile for arccfd codes                      #
#                                                                             #
###############################################################################


all: library documentation


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


clean: cleanCode cleanDoc cleanLib
	rm -f *~ */*~


cleanCode:
	@echo ""
	@echo "cleaning code . . ."
	@echo ""
	rm -rf include
	cd src; make clean


cleanDoc:
	@echo ""
	@echo "cleaning documentation . . ."
	@echo ""
	cd doc; rm -rf html latex *~


cleanLib:
	@echo ""
	@echo "cleaning library . . ."
	@echo ""
	cd lib; rm -rf *.a *~
