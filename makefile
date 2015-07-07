if [[-z "$CXX"]]; then CXX='c++'
CFLAGS=-c -std=c++11
LIBS=
VPATH=./src
OBJECTS_NEUTRAL = ClusterData.o Utils.o ClusterCalculations.o mainNeutral.o
OBJECTS_ACTUAL = ClusterData.o Utils.o ClusterCalculations.o mainActual.o

LL=latex
PDF=dvipdfm
LFLAGS=-interaction=nonstopmode -shell-escape
LFILE = aaaread-this

# ------ C++ make

all: mainNeutral mainActual

neutral: mainNeutral

actual: mainActual

mainNeutral: $(OBJECTS_NEUTRAL)
	$(CXX) $(OBJECTS_NEUTRAL) -o ClustersNeutral $(LIBS) 

mainActual: $(OBJECTS_ACTUAL)
	$(CXX) $(OBJECTS_ACTUAL) -o ClustersActual $(LIBS) 

%.o: %.c++
	$(CXX) $(CFLAGS) $<

# ------ latex make

#latex: latex1 latex1 dvi clean
latex: latex1 latex1 pdf clean

latex1: 
	$(LL) $(LFILE).tex $(LFLAGS) 

pdf:
	$(PDF) $(LFILE).dvi

dvi:
	dvips -P pdf -q $(LFILE).dvi
	ps2pdf $(LFILE).ps

clean:
	rm -f *.o *.dvi *.blg *.maf *.mtc *.mtc0 *.ps *.soc

