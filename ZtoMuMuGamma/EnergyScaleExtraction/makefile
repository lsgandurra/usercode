# définition des cibles particulières
.PHONY: clean, mrproper
 
# désactivation des règles implicites
.SUFFIXES:
 
# définition des variables
CC = g++
#CFLAGS = -W -Wall -v -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags`
#CFLAGS = -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags` 
CFLAGS = `root-config --cflags --libs` -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -W 
  
# all
all: SFits.o fitFunctions.o functions.o
	$(CC) SFits.o fitFunctions.o functions.o -o SFits.exe $(CFLAGS)
 
SFits.o: SFits.cpp fitFunctions.h functions.h
	$(CC) -c SFits.cpp -o SFits.o $(CFLAGS)

fitFunctions.o : fitFunctions.cpp
	$(CC) -c fitFunctions.cpp -o fitFunctions.o $(CFLAGS)

functions.o : functions.cpp
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
 
# clean
clean:
	rm -rf *.bak rm -rf *.o
 
# mrproper
mrproper: clean
	rm -rf Programme





