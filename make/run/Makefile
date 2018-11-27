include ../*.mk

all: $(MAIN)
.PHONY: all

## check and recompile all *.o files upon which MAIN depends
$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

## link *.o files to create MAIN
$(MAIN): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

## prevent creation of CLEAN file
.PHONY: cleanall cleanobj clean~ cleanprog

## clean directories of *.o and ~ files
cleanall: cleanobj clean~ cleanprog

cleanobj:
	rm -f $(ODIR)/*.o

clean~:
	rm -f *~

cleanprog:
	rm -f $(MAIN)