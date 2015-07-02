if [[-z "$CXX"]]; then CXX='c++'
CFLAGS=-c -std=c++11
LIBS=
VPATH=./src
OBJECTS_NEUTRAL = ClusterData.o Utils.o ClusterCalculations.o mainNeutral.o
OBJECTS_ACTUAL = ClusterData.o Utils.o ClusterCalculations.o mainActual.o

all: mainNeutral mainActual

neutral: mainNeutral

actual: mainActual

mainNeutral: $(OBJECTS_NEUTRAL)
	$(CXX) $(OBJECTS_NEUTRAL) -o ClustersNeutral $(LIBS) 

mainActual: $(OBJECTS_ACTUAL)
	$(CXX) $(OBJECTS_ACTUAL) -o ClustersActual $(LIBS) 

%.o: %.c++
	$(CXX) $(CFLAGS) $<

clean:
	rm -f *.o 

