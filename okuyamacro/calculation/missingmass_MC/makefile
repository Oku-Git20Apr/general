TARGET1=missingmass_K
OBJS1=missingmass_K.o
TARGET2=missingmass_Mg
OBJS2=missingmass_Mg.o
TARGET3=missingmass_B_e05
OBJS3=missingmass_B_e05.o
TARGET2=missingmass_Mg_highq
OBJS2=missingmass_Mg_highq.o
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

#.PHONY: all
all: $(TARGET1)\
	$(TARGET2)\
	$(TARGET3)\
	$(TARGET4)

#test2: test2.o
#	g++ -o $@ test2.o $(CXXLIBS)
#
#.cc.o:
#	g++ -c $(CXXFLAGS) $<

${TARGET1}: ${OBJS1}
	g++ -o $(TARGET1) ${OBJS1} $(CXXLIBS) $(CXXFLAGS)
${TARGET2}: ${OBJS2}
	g++ -o $(TARGET2) ${OBJS2} $(CXXLIBS) $(CXXFLAGS)
${TARGET3}: ${OBJS3}
	g++ -o $(TARGET3) ${OBJS3} $(CXXLIBS) $(CXXFLAGS)
${TARGET4}: ${OBJS4}
	g++ -o $(TARGET4) ${OBJS4} $(CXXLIBS) $(CXXFLAGS)

#.PHONY: clean
clean:
	rm -f *.o core ${TARGET1} ${TARGET2} ${TARGET3} ${TARGET4}
