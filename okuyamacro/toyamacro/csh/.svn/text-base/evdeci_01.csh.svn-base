#!/bin/csh -f

  if ( $#argv < 2 ) then
    echo ""
    echo "  Usage: ./evdeci.csh [period] [run number]"
    echo ""
    exit
  endif

  set period = ${1}
  set run    = ${2}

  set user_name = `whoami`

  set base_dir  = `pwd`/../../../


  set log_dir  = ${base_dir}job/log/EvnumDeci/${run}/
  
  set host_num = `echo ${HOST} | sed 's/farm//g'`
  set data_dir = /data/${host_num}c
  if (! -e ${data_dir} ) then
  set data_dir = /data/${host_num}b
  endif
  if (! -e ${data_dir} ) then
  set data_dir = /data/${host_num}d
  endif

  set job_tmp_dir  = ${data_dir}/.condor_job_tmp/${user_name}/EvnumDeci_${period}_${run}

  if (! -e ${log_dir} ) then
    mkdir -p ${log_dir}
  endif

  if (! -e ${job_tmp_dir} ) then
    mkdir -p ${job_tmp_dir}
  endif

  set log_file    =     ${log_dir}/EvnumDeci_${period}_${run}.log
  set conf_file   = ${job_tmp_dir}/EvnumDeci_${period}_${run}.conf
  set output_file = ${job_tmp_dir}/dummy.root


  set exec_file = ${base_dir}bin/EvnumDecimate

cat << EOF >>  ${conf_file}
CMAP:	param/countermap_from10068.param
WMAP:	param/wiremap.param
FIELD:	ToscaField/680MagDef.map
NOMAG:  0
POL:    -1
DCGEO:	param/dcgeo.param
DCTDC:	param/dctdc/10246.param
DCDRFT:	param/dcdrift.param_calib_step_03
HDGEO:	param/hodogeo.param
HDCOR:	param/hodotune_de/${run}.param
HDPHC:  param/hodophc_10246_itr0.param
HDPHC:  param/hodophc_10246_itr1.param
HDPHC:  param/hodophc_10246_itr2.param
HDPHC:  param/hodophc_10246_itr3.param
TAGE:	param/tagger.param
TAGT:	param/taggertime.param
TAGA:   param/tagb.param
DAQ3T:  param/TaggerTime/10246.param
DAQ3E:  param/TaggerInfo.param
SCA:	param/scalerdef.param
TARGP:	0.0  0.0  0.0
GTAG:	-2. 2. -2. 2.
EOF

cat ${base_dir}/macros/toyamacro/conf/run${run}plus_minus_event.dat_uniq_sort >>  ${conf_file}

#cat ${conf_file}
  cp  ${base_dir}k0daq0/r00${run}.dat  ${job_tmp_dir}/k0daq0_r00${run}.dat
  cp  ${base_dir}k0daq1/r00${run}.dat  ${job_tmp_dir}/k0daq1_r00${run}.dat
  cp  ${base_dir}k0daq2/r00${run}.dat  ${job_tmp_dir}/k0daq2_r00${run}.dat
  cp  ${base_dir}k0daq3/r00${run}.dat  ${job_tmp_dir}/k0daq3_r00${run}.dat
  chmod 644                                   ${job_tmp_dir}/*.dat


  date >& ! ${log_file}
  # echo $LD_LIBRARY_PATH >& ! ${log_file}

  cd ${base_dir}
  ${exec_file} ${conf_file} ${output_file} ${job_tmp_dir}/k0daq0_r00${run}.dat ${job_tmp_dir}/k0daq1_r00${run}.dat ${job_tmp_dir}/k0daq2_r00${run}.dat ${job_tmp_dir}/k0daq3_r00${run}.dat \
                                           ${job_tmp_dir}/k0daq0_${run}_deci.dat ${job_tmp_dir}/k0daq1_${run}_deci.dat ${job_tmp_dir}/k0daq2_${run}_deci.dat ${job_tmp_dir}/k0daq3_${run}_deci.dat >>& ${log_file}

  echo mv ${job_tmp_dir}/k0daq0_${run}_deci.dat  ${base_dir}k0daq0_pm/ >>& ${log_file}
       mv ${job_tmp_dir}/k0daq0_${run}_deci.dat  ${base_dir}k0daq0_pm/ >>& ${log_file}
  echo mv ${job_tmp_dir}/k0daq1_${run}_deci.dat  ${base_dir}k0daq1_pm/ >>& ${log_file}
       mv ${job_tmp_dir}/k0daq1_${run}_deci.dat  ${base_dir}k0daq1_pm/ >>& ${log_file}
  echo mv ${job_tmp_dir}/k0daq2_${run}_deci.dat  ${base_dir}k0daq2_pm/ >>& ${log_file}
       mv ${job_tmp_dir}/k0daq2_${run}_deci.dat  ${base_dir}k0daq2_pm/ >>& ${log_file}
  echo mv ${job_tmp_dir}/k0daq3_${run}_deci.dat  ${base_dir}k0daq3_pm/ >>& ${log_file}
       mv ${job_tmp_dir}/k0daq3_${run}_deci.dat  ${base_dir}k0daq3_pm/ >>& ${log_file}

  echo rm -f  ${job_tmp_dir}/${output_file} >>& ${log_file}
       rm -f  ${job_tmp_dir}/${output_file} >>& ${log_file}
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

