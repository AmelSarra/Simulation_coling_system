CFLAGS = -O
CC = g++

Projet : Projet.o Ailette.o functions_stationary.o functions_instationary.o functions_3d.o
	$(CC) $(CFLAGS) -o Projet Projet.o Ailette.o functions_stationary.o functions_instationary.o functions_3d.o
	make clean

Projet.o:
	$(CC) $(CFLAGS) -c src/Projet.cpp

Ailette.o:
	$(CC) $(CFLAGS) -c src/classes/Ailette.cpp

functions_stationary.o:
	$(CC) $(CFLAGS) -c src/functions/functions_stationary.cpp

functions_instationary.o:
	$(CC) $(CFLAGS) -c src/functions/functions_instationary.cpp

functions_3d.o:
	$(CC) $(CFLAGS) -c src/functions/functions_3d.cpp

clean:
	rm -f core *.o
