#!/bin/csh -f

  if ( $#argv < 4 ) then
    echo ""
    echo "  Usage: ./rkvertex.csh [period] [run number] [begin event] [end event]"
    echo ""
    exit
  endif

  set period = ${1}
  set run    = ${2}
  set bgnevt = ${3}
  set endevt = ${4}

  set user_name = `whoami`

  set base_dir  = `pwd`/../../../


  set log_dir  = ${base_dir}job/log/RKVertex/${run}/
  set root_dir = ${base_dir}job/root/RKVertex/${run}/
  
  set host_num = `echo ${HOST} | sed 's/farm//g'`
  #echo ${host_num}
  set data_dir = /data/${host_num}c
  if (! -e ${data_dir} ) then
  set data_dir = /data/${host_num}b
  else if (${host_num}>12 && ${host_num}<14) then
  set data_dir = /data/${host_num}b
  endif
  if (! -e ${data_dir} ) then
  set data_dir = /data/${host_num}d
  endif

  set job_tmp_dir  = ${data_dir}/.condor_job_tmp/${user_name}/RKTrack_${period}_${run}_${bgnevt}_${endevt}
  ##set job_tmp_dir  = /data/${host_num}c/.condor_job_tmp/${user_name}/RKVertex_${period}_${run}_${bgnevt}_${endevt}

  if (! -e ${log_dir} ) then
    mkdir -p ${log_dir}
  endif

  if (! -e ${root_dir} ) then
    mkdir -p ${root_dir}
  endif

  if (! -e ${job_tmp_dir} ) then
    mkdir -p ${job_tmp_dir}
  endif

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

  set log_file    =     ${log_dir}/RKVertex_${period}_${run}${seven_digit_bgnevt}${seven_digit_endevt}.log
  set conf_file   = ${job_tmp_dir}/RKVertex_${period}_${run}${seven_digit_bgnevt}${seven_digit_endevt}.conf
  set output_file = ${job_tmp_dir}/RKVertex_${period}_${run}${seven_digit_bgnevt}${seven_digit_endevt}.root


  set exec_file = ${base_dir}bin/RKVertexAnalysis

# cat ${base_dir}/conf/${period}/${run}.conf >!  ${conf_file}
#cat ${base_dir}/conf/${period}/analyzer.conf >!  ${conf_file}
#cat ${base_dir}conf/10246.conf >!  ${conf_file}
cat << EOF >>!  ${conf_file}
CMAP:	param/countermap_from10068.param
WMAP:	param/wiremap.param
FIELD:	ToscaField/680MagDef.map
NOMAG:  0
POL:    -1
DCGEO:	param/dcgeo.param
DCTDC:	param/dctdc/10246.param
DCDRFT:	param/dcdrift.param_calib_step_03
HDGEO:	param/hodogeo.param
HDCOR:	param/hodo_combine/${run}.param
HDPHC:  param/hodophc_10246_itr0.param
HDPHC:  param/hodophc_10246_itr1.param
HDPHC:  param/hodophc_10246_itr2.param
#HDPHC:  param/hodophc_10246_itr3.param
TAGE:	param/tagger.param
TAGT:	param/taggertime.param
TAGA:   param/tagb.param
DAQ3T:  param/TaggerTime/10246_01.param
#DAQ3T:  param/TaggerTime/${run}.param
DAQ3E:  param/TaggerInfo.param
SCA:	param/scalerdef.param
TARGP:	0.0  0.0  0.0
GTAG:	-2. 2. -2. 2.
BEGEV:  ${bgnevt}
ENDEV:  ${endevt}
EOF


  cp  ${base_dir}k0daq0/r00${run}.dat  ${job_tmp_dir}/k0daq0_r00${run}.dat
  cp  ${base_dir}k0daq1/r00${run}.dat  ${job_tmp_dir}/k0daq1_r00${run}.dat
  cp  ${base_dir}k0daq2/r00${run}.dat  ${job_tmp_dir}/k0daq2_r00${run}.dat
  cp  ${base_dir}k0daq3/r00${run}.dat  ${job_tmp_dir}/k0daq3_r00${run}.dat
  chmod 644                                   ${job_tmp_dir}/*.dat


  date >& ! ${log_file}
  # echo $LD_LIBRARY_PATH >& ! ${log_file}

  cd ${base_dir}
  ${exec_file} ${conf_file} ${output_file} ${job_tmp_dir}/k0daq0_r00${run}.dat ${job_tmp_dir}/k0daq1_r00${run}.dat ${job_tmp_dir}/k0daq2_r00${run}.dat ${job_tmp_dir}/k0daq3_r00${run}.dat >>& ${log_file}


  echo mv ${output_file}  ${root_dir} >>& ${log_file}
       mv ${output_file}  ${root_dir} >>& ${log_file}

  echo rm -f  ${job_tmp_dir}/k0daq0_r00${run}.dat >>& ${log_file}
       rm -f  ${job_tmp_dir}/k0daq0_r00${run}.dat >>& ${log_file}
  echo rm -f  ${job_tmp_dir}/k0daq1_r00${run}.dat >>& ${log_file}
       rm -f  ${job_tmp_dir}/k0daq1_r00${run}.dat >>& ${log_file}
  echo rm -f  ${job_tmp_dir}/k0daq2_r00${run}.dat >>& ${log_file}
       rm -f  ${job_tmp_dir}/k0daq2_r00${run}.dat >>& ${log_file}
  echo rm -f  ${job_tmp_dir}/k0daq3_r00${run}.dat >>& ${log_file}
       rm -f  ${job_tmp_dir}/k0daq3_r00${run}.dat >>& ${log_file}

  echo rm -rf ${job_tmp_dir}/ >>& ${log_file}
       rm -rf ${job_tmp_dir}/

  
  date >>&  ${log_file}

