LIB        = -L. 
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = MD.exe
CXX        = g++

${EXEC}: MD.c
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} MD.c -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<
