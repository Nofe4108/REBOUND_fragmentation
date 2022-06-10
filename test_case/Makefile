export OPENGL=1

REB_DIR=/Users/nferich/GitHub/rebound

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I /Users/nferich/GitHub/rebound/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C /Users/nferich/GitHub/rebound/src/
	@-rm -f librebound.so
	@ln -s /Users/nferich/GitHub/rebound/src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C /Users/nferich/GitHub/rebound/src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
