CC=g++
CFLAGS=-c
LIBS=-lboost_system
VPATH=./src
OBJECTS1 = NeutralClusters.o Utils.o InOut.o DataProcessing.o Calculations.o
OBJECTS2 = GetR2Mats.o Utils.o InOut.o DataProcessing.o Calculations.o
OBJECTS3 = Clusters.o Utils.o InOut.o DataProcessing.o Calculations.o

%.o: %.c++
	$(CC) $(CFLAGS) $<

all: ncl

ncl: $(OBJECTS1)
	$(CC) $(OBJECTS1) $(LIBS) -o NeutralClusters

getr2: GetR2Mats

GetR2Mats: $(OBJECTS2)
	$(CC) $(OBJECTS2) $(LIBS) -o getr2

clusters: Clusters

Clusters: $(OBJECTS3)
	$(CC) $(OBJECTS3) $(LIBS) -o Clusters

clean:
	rm -f *.o 
