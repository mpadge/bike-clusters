CC=g++
CFLAGS=-c
LIBS=
VPATH=./src
OBJECTS = NeutralClusters.o Utils.o InOut.o DataProcessing.o Calculations.o

all: ncl

ncl: $(OBJECTS)
	$(CC) $(OBJECTS) -o $(LIBS) NeutralClusters

$(BUILDDIR)/%.o: $(BUILDDIR)/%.c++
	$(CC) $(CFLAGS) $@ $< -o

clean:
	rm -f *.o 

# ---------------------------------------------------

OBJECTS2 = GetR2Mats.o Utils.o InOut.o DataProcessing.o Calculations.o

getr2: GetR2Mats

GetR2Mats: $(OBJECTS2)
	$(CC) $(OBJECTS2) -o $(LIBS) getr2

%.o: %.c++
	$(CC) $(CFLAGS) $<

# ---------------------------------------------------

OBJECTS3 = Clusters.o Utils.o InOut.o DataProcessing.o Calculations.o

clusters: Clusters

Clusters: $(OBJECTS3)
	$(CC) $(OBJECTS3) -o $(LIBS) Clusters

%.o: %.c++
	$(CC) $(CFLAGS) $<

