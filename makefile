CC=g++ 
OBJ=	main.o\
CC=g++
#FFTW_INC_INCLUDE = /Users/JZ/Desktop/ComMathPhys/Comp/fftw-3.3.3/include
FFTW_INC_LIB  = /usr/local/include
OBJ=    main.o\

PSTD: $(OBJ)
	$(CC) -o cpstd  $(OBJ) -L $(FFTW_INC_LIB) -lfftw3 -lm
#       $(C) -o cpstd  $(OBJ)  -lfftw3 -lm
main.o: main.cpp
	$(CC) -c  main.cpp
#initialization.o:initialization.cpp
#	$(CC) -c  initialization.cpp
#simulation.o:simulation.cpp
#	$(CC) -c  simulation.cpp
#transform.o:transform.cpp
#	$(CC) -c  transform.cpp -I $(FFTW_INC_LIB)
#source.o:source.cpp
#	$(CC) -c  source.cpp
#memory.o:memory.cpp
#	$(CC) -c  memory.cpp
clean::
	rm -f *.o  cpstd *.txt

