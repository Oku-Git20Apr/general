#!/bin/csh
  set period = 2017_04

  set user_name = `whoami`
  set base_dir  = `pwd`/../../../
  set incevt    = 20000
  date >> log/rkvertex01.log

#foreach run (10246)
#foreach run (10246 10247 10248 10249 10250 10251 10252 10253 10254 10255 \
#             10256 10257 10258 10259 10261 10262 10268 10269 10270 10271 \
#             10272 10273 10274 10275 10276 10277 10278 10280 10282 10283 \
#             10284 10286 10287 10288 10289 10290 10291 10292 \
#)
#foreach run (10226 )
#             10225 10226 \
#             10227 10228 10229 10230 10232 10233 10234 10236 10237 10238 \
#             10239 \
#foreach run ( \
#             10226 \
#             10238 \
#             10239 10240 10241 10242 10243 10245 \
#)
#             10293 10294 \
#             10296 10297 10298 10304 10305 10307 10308 10309 10309 10310 \
#             10311 10312 10313 10314 10315 10316 10317 10318 10320 10321 \
#             10322 10323 10324 10324 10325 10326 10327 10342 10345 10346 \
#             10347 10348 10350 )
# foreach run ( `cut -c-95 ${base_dir}/data/run_summary_${period}.txt | grep -i normal | tr 'n' ' ' | awk '{print $2}'` )


#foreach run (10246 10257 10270 10282 10294 10310 10321) ##standard runs

##group1##st10226 
##group2##st10238 
#foreach run (10226 10227 10228 \
#             10229 10230 10232 10233 \
#             10234 10236 10237 10238 10239 10240 10241 \
#             10242 10243 10245 )

##group2##st10238 
#foreach run (10239 10240 10241 \
#             10242 10243 10245  )

##group3##st10246 
#foreach run (10246 10247 10248 10249 \
#             10250 10251 10252 10253 \
#             10254 10255 10256 \
#             10258 10259 10261  \
#             10262 10268 10269)
##group4##st10257  
#foreach run (10258 10261  \
#             10262 10268 10269 )

##group5##st10270 
#foreach run (10270 10271 10272 10273  \
#             10274 10275 10276 10277  \
#             10278 10280  )
#foreach run (10257 10270)

##group6##st10282 
foreach run (10282 10283 10284 10286  \
             10287 10288 10289 10290  \
             10291 10292 10293 10294 10296 10297 10298  \
             10304 10305 10307 10308  \
             10309 10310\
             10311 10312 10313  \
             10314 10315 10316 10317  \
             10318 10320 )

##group7##st10294 
#foreach run (10296 10297 10298  \
#             10304 10305 10307 10308  \
#             10309  \
#             10311 10312 10313  \
#             10314 10315 10316 10317  \
#             10318 10320 \
#             )

##group8##st10310 
#foreach run (10304 10305 10307 10308  \
#             10309  \
#             10311 10312 10313  \
#             10314 10315 10316 10317  \
#             10318 10320  )
##group9##st10321 
#foreach run (10321 10322 10323 10324  \
#             10325 10326 10342 10345  \
#             10346 10347 10348 10350  )
#foreach run (10246 10248 )

####OHH param tune
#foreach run (10247 10248 10249 10250 10251 10252)
foreach run (10259)

    set bgnevt = 1
    @ endevt = ${bgnevt} + ${incevt}
    @ endevt = ${endevt} - 1

    set root_dir = ${base_dir}job/root/RKVertex/${run}

    set nevent = `cut -c-95 ${base_dir}data/run_summary_${period}.txt |grep run${run} | awk '{print $3}'`
    while ( ${bgnevt} < ${nevent} )

echo ${run} ${bgnevt} ${nevent} ${nevent}
      if ( ${bgnevt} < 10 ) then
        set seven_digit_bgnevt = _000000${bgnevt}
      else if ( ${bgnevt} < 100 ) then
        set seven_digit_bgnevt = _00000${bgnevt}
      else if ( ${bgnevt} < 1000 ) then
        set seven_digit_bgnevt = _0000${bgnevt}
      else if ( ${bgnevt} < 10000 ) then
        set seven_digit_bgnevt = _000${bgnevt}
      else if ( ${bgnevt} < 100000 ) then
        set seven_digit_bgnevt = _00${bgnevt}
      else if ( ${bgnevt} < 1000000 ) then
        set seven_digit_bgnevt = _0${bgnevt}
      else
        set seven_digit_bgnevt = _${bgnevt}
      endif

      if ( ${endevt} < 10 ) then
        set seven_digit_endevt = _000000${endevt}
      else if ( ${endevt} < 100 ) then
        set seven_digit_endevt = _00000${endevt}
      else if ( ${endevt} < 1000 ) then
        set seven_digit_endevt = _0000${endevt}
      else if ( ${endevt} < 10000 ) then
        set seven_digit_endevt = _000${endevt}
      else if ( ${endevt} < 100000 ) then
        set seven_digit_endevt = _00${endevt}
      else if ( ${endevt} < 1000000 ) then
        set seven_digit_endevt = _0${endevt}
      else
        set seven_digit_endevt = _${endevt}
      endif

      set output_file = ${root_dir}/RKVertex_${period}_${run}${seven_digit_bgnevt}${seven_digit_endevt}.root

      if (-e ${output_file} ) then
        echo ${output_file} exist.
        set file_size = `ls -l ${output_file} | awk '{print $5}'`
        echo "size of root file is ${file_size}"
          if(${file_size}>10000) then
            echo "file size is reasonable (more than 10k)"
          else
            echo "file size is too small (less than 10k)"
          endif
      else
        echo No ${output_file}.
        set file_size = 0
      endif

  cat << EOF >!  .rkvertex_01.input
  Universe     = vanilla
  Initialdir   = ${base_dir}macros/toyamacro/job
  Executable   = ${base_dir}macros/toyamacro/csh/rkvertex_01.csh
  Arguments    = ${period} ${run} ${bgnevt} ${endevt}
  +MachineList0  = Machine != "farm1.npex.prv"
  +MachineList1  = Machine != "farm2.npex.prv"
  +MachineList2 = Machine != "farm4.npex.prv"
  +MachineList3 = Machine != "farm5.npex.prv"
  +MachineList4 = Machine != "farm6.npex.prv"
  +MachineList5 = Machine != "farm12.npex.prv"
  +MachineList6 = Machine != "farm13.npex.prv"
  +MachineList7 = Machine != "farm19.npex.prv"
  +MachineList8 = Machine != "farm30.npex.prv"
  +MachineList9 = Machine != "farm31.npex.prv"
  +MachineList  = MachineList0 && MachineList1 && MachineList2 && MachineList3 && MachineList4 && MachineList5 && MachineList6 && MachineList7 && MachineList8 && MachineList9
  Requirements = MachineList && ( Arch == "X86_64" ) && (OpSys == "LINUX") && ( Memory > 500 )
  getenv       = True
  notification = Error
  notify_user  = ${user_name}@lambda.phys.tohoku.ac.jp
  Queue
EOF


set n=1
while(${n} <= 1)
  set tmp_job_num_run  = `condor_status -submitter |grep ${user_name} | head -n 1 | awk '{print $3}'`
  set tmp_job_num_idle = `condor_status -submitter |grep ${user_name} | head -n 1 | awk '{print $4}'`
  set job_num_run = 0
  set job_num_idle = 0
  set job_num_total = 0
  @ job_num_run   = ${tmp_job_num_run} + ${job_num_run}
  @ job_num_idle  = ${tmp_job_num_idle} + ${job_num_idle}
  echo ${job_num_run} ${job_num_idle}


  @ job_num_total = ${job_num_run} + ${job_num_idle}
  set job_num = ${job_num_idle}
  #set job_num = `condor_q | grep -c ${user_name} | awk '{print $1}'`
  echo ${job_num}   `date`
  
  if(${job_num}<100)then
     @ n = ${n} + 1
  else
    echo "sleeping..."
    sleep 20m
  endif
end
  
if ( (-e ${output_file} && ${file_size}<10000) || ! -e ${output_file} ) then
          condor_submit .rkvertex_01.input
  echo "submit ${run}  ${bgnevt} ${endevt}" >> log/rkvertex01.log
  
else
  echo "root file is ok. ${run}"
endif

 
      @ bgnevt = ${bgnevt} + ${incevt}
      @ endevt = ${endevt} + ${incevt}

    end


 end
