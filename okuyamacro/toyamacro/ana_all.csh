#!/bin/csh

set dtime = "4ns"
foreach run (   \
#             10226_10245 \
             10226_10238 \
             10246_10256 \
             10257_10269 \
             10270_10280 \
             10281_10292 \
             10294_10309 \
             10310_10320 \
             10321_10350 \
            )

./bin/RKVertex_ana -f ./list/${run}_t0.dat -w root/RKV_ana_acci/ana${run}_${dtime}.root -p pdf/vertex/ana${run}_${dtime}.pdf 
#./bin/RKVertex_ana -f ./list/${run}_t0.dat -w root/ana${run}.root -p pdf/vertex/ana${run}.pdf 
#./bin/RKVertex_ana -f ./list/${run}_pm.dat -w root/RKV_ana_pm/ana${run}.root -p pdf/vertex/ana${run}.pdf &

end

#sleep 5m
#./hadd_ana.csh

