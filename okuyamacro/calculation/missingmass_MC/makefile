TARGET1=missingmass_K
OBJS1=missingmass_K.o
TARGET2=missingmass_Mg
OBJS2=missingmass_Mg.o
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

#.PHONY: all
all: $(TARGET1)\
	$(TARGET2)

#test2: test2.o
#	g++ -o $@ test2.o $(CXXLIBS)
#
#.cc.o:
#	g++ -c $(CXXFLAGS) $<

${TARGET1}: ${OBJS1}
	g++ -o $(TARGET1) ${OBJS1} $(CXXLIBS) $(CXXFLAGS)
${TARGET2}: ${OBJS2}
	g++ -o $(TARGET2) ${OBJS2} $(CXXLIBS) $(CXXFLAGS)

#.PHONY: clean
clean:
	rm -f *.o core ${TARGET1} ${TARGET2}
