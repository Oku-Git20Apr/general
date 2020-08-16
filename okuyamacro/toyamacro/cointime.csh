#!/bin/csh

if ( $#argv < 1 ) then
  echo ""
  echo "  Usage: ./cointime.csh [5 digit run number]"
  echo ""
  exit
endif


set run = ${1}

###./bin/Cointime_TDL_TagB -f ../../root/all/${1}.root -w root/${1}_coin.root
./bin/tdl_ana1 -f ../../root/all/${1}.root -w root/${1}_coin.root

###TDL raw TDC data
root -b << EOF
.L draw_TDLRaw.C
draw_TDLRaw("${run}_coin.root");
EOF

##coin time TagB Sum vs TDL 2-9
root -b << EOF
.L draw_TagSum.C
draw_TagSum("${run}_coin.root");
EOF

###coin time TagB25 vs TDL 2-9
#root -b << EOF
#.L draw_TagB25.C
#draw_TagB25("${run}_coin.root");
#EOF
#
###coin time TagB 1-30 vs TDL 2
#root -b << EOF
#.L draw_TDL2.C
#draw_TDL2("${run}_coin.root");
#EOF
#
###coin time TagB 25 vs TDL R2 evaluate signal to noise ratio
#root -b << EOF
#.L TagB25_TDLR2SN.C
#TagB25_TDLR2SN("${run}_coin.root");
#EOF

##TDL ToF resolution
#root -b << EOF
#.L draw_TekitoToF.C
#draw_TekitoToF("${run}_coin.root");
#EOF


#./bin/TDLMonitor -f ../../root/all/${1}.root -w pdf/${1}_tdlmoni.pdf
