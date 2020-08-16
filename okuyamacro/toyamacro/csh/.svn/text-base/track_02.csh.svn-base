#!/bin/csh -f

  if ( $#argv < 4 ) then
    echo ""
    echo "  Usage: ./track_02.csh [root_dir] [root_file_base] [i_first] [i_end]"
    echo "[i_first] : first file number to combine"
    echo "[i_end] : last file number  to combine"
    exit
  endif

  set root_dir       = ${1}
  set root_file_base = ${2}
  set i_first        = ${3}
  set i_end          = ${4}

  @ ntail = ${i_end} - ${i_first}
  @ ntail = ${ntail} + 1

  @ nhead = ${i_end} 

  if ( ${i_first} < 10 ) then
    set isub = _00${i_first}
  else if ( ${i_first} < 100 ) then
    set isub = _0${i_first}
  else
    set isub = _${i_first}
  endif


  cd ${root_dir}

  if (! -e  ${root_file_base}${isub}.root ) then
    hadd ${root_file_base}${isub}.root `ls | grep ${root_file_base}_ | head -${nhead} | tail -${ntail}`
  endif



