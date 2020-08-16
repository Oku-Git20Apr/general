CC=gcc
CXX=g++
#CC=icc
#CXX=icpc
#CFLAGS  = -parallel -par-report3 -par-threshold0 -O3

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

TARGET1=     tdl_ana1
OBJS1=       tdl_ana1.o
TARGET2=     cre_vmethre
OBJS2=       cre_vmethre.o
TARGET3=     ped_supp
OBJS3=       ped_supp.o
TARGET4=     chain
OBJS4=       chain.o
TARGET5=     gain_check
OBJS5=       gain_check.o
TARGET6=     TDLreso_V1290
OBJS6=       TDLreso_V1290.o
TARGET7=     Cointime_TDL_TagB
OBJS7=       Cointime_TDL_TagB.o
TARGET8=     TDLMonitor
OBJS8=       TDLMonitor.o
TARGET9=     tdl_ana2
OBJS9=       tdl_ana2.o Tree.o
TARGET10=    track_ana
OBJS10=      track_ana.o Track_tr.o
TARGET11=    InvariantMass
OBJS11=      InvariantMass.o
TARGET12=    proton_select_wotrack
OBJS12=      proton_select_wotrack.o Tree.o
TARGET13=    RKVertex_ana
OBJS13=      RKVertex_ana.o RKVertex_tr.o Setting.o ToF_ParamMan.o
TARGET14=    RKVertex_remakeroot
OBJS14=      RKVertex_remakeroot.o RKVertex_tr.o RKV_ParamMan.o
TARGET15=    RKVertex_mkparam
OBJS15=      RKVertex_mkparam.o RKVertex_tr.o Setting.o RKV_ParamMan.o
TARGET16=    RKVertex_ppi_extract
OBJS16=      RKVertex_ppi_extract.o RKVertex_tr.o RKV_ParamMan.o
TARGET17=     RF_ana
OBJS17=       RF_ana.o Tree.o Setting.o
TARGET18=     combine_param
OBJS18=       combine_param.o HodoParamMan.o RKV_ParamMan.o TagParamMan.o
TARGET19=     phc_check
OBJS19=       phc_check.o Tree.o Setting.o
TARGET20=     RKTrack_geocheck
OBJS20=       RKTrack_geocheck.o RKTrack_tr.o Setting.o HodoGeomMan.o
TARGET21=    RKVertex_modparam
OBJS21=      RKVertex_modparam.o RKV_ParamMan.o
TARGET22=    date
OBJS22=      date.o Setting.o
TARGET23=    RKVertex_selectevent
OBJS23=      RKVertex_selectevent.o RKVertex_tr.o
TARGET24=    RKVertex_mkparam_ft_corr
OBJS24=      RKVertex_mkparam_ft_corr.o RKVertex_tr.o Setting.o ToF_ParamMan.o
TARGET25=    draw_tdlparam
OBJS25=      draw_tdlparam.o RKV_ParamMan.o Setting.o
TARGET26=    check_coin
OBJS26=      check_coin.o HodoParamMan.o RKV_ParamMan.o TagParamMan.o Setting.o
TARGET27=    checkparam
OBJS27=      checkparam.o HodoParamMan.o
TARGET28=    minuit
OBJS28=      minuit.o MathTools.o

all: $(TARGET1) \
     $(TARGET2) \
     $(TARGET3) \
     $(TARGET4) \
     $(TARGET5) \
     $(TARGET6) \
     $(TARGET7) \
     $(TARGET8) \
     $(TARGET9) \
     $(TARGET10) \
     $(TARGET11) \
     $(TARGET12) \
     $(TARGET13) \
     $(TARGET14) \
     $(TARGET15) \
     $(TARGET16) \
     $(TARGET17) \
     $(TARGET18) \
     $(TARGET19) \
     $(TARGET20) \
     $(TARGET21) \
     $(TARGET22) \
     $(TARGET23) \
     $(TARGET24) \
     $(TARGET25) \
     $(TARGET26) \
     $(TARGET27) \
     $(TARGET28) \

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET5): $(patsubst %,$(LIBDIR)/%,$(OBJS5))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET6): $(patsubst %,$(LIBDIR)/%,$(OBJS6))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET7): $(patsubst %,$(LIBDIR)/%,$(OBJS7))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET8): $(patsubst %,$(LIBDIR)/%,$(OBJS8))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET9): $(patsubst %,$(LIBDIR)/%,$(OBJS9))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET10): $(patsubst %,$(LIBDIR)/%,$(OBJS10))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET11): $(patsubst %,$(LIBDIR)/%,$(OBJS11))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET12): $(patsubst %,$(LIBDIR)/%,$(OBJS12))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET13): $(patsubst %,$(LIBDIR)/%,$(OBJS13))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET14): $(patsubst %,$(LIBDIR)/%,$(OBJS14))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET15): $(patsubst %,$(LIBDIR)/%,$(OBJS15))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET16): $(patsubst %,$(LIBDIR)/%,$(OBJS16))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET17): $(patsubst %,$(LIBDIR)/%,$(OBJS17))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET18): $(patsubst %,$(LIBDIR)/%,$(OBJS18))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET19): $(patsubst %,$(LIBDIR)/%,$(OBJS19))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET20): $(patsubst %,$(LIBDIR)/%,$(OBJS20))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET21): $(patsubst %,$(LIBDIR)/%,$(OBJS21))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET22): $(patsubst %,$(LIBDIR)/%,$(OBJS22))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET23): $(patsubst %,$(LIBDIR)/%,$(OBJS23))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET24): $(patsubst %,$(LIBDIR)/%,$(OBJS24))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET25): $(patsubst %,$(LIBDIR)/%,$(OBJS25))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET26): $(patsubst %,$(LIBDIR)/%,$(OBJS26))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET27): $(patsubst %,$(LIBDIR)/%,$(OBJS27))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET28): $(patsubst %,$(LIBDIR)/%,$(OBJS28))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*

