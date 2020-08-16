#!/bin/csh -f

  if ( $#argv < 2 ) then
    echo ""
    echo "  Usage: ./RKVertex_05.csh [root_dir] [root_file_base]"
    echo ""
    exit
  endif

  set root_dir       = ${1}
  set root_file_base = ${2}

  cd ${root_dir}

  if (-e ${root_file_base}.root ) then
    rm -f  `ls ${root_file_base}_*.root`
  endif

