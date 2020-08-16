#!/bin/csh -f

  if ( $#argv < 2 ) then
    echo ""
    echo "  Usage: ./RKVertex_04.csh [root_dir] [root_file_base]"
    echo ""
    exit
  endif

  set root_dir       = ${1}
  set root_file_base = ${2}

  cd ${root_dir}

  hadd -f ${root_file_base}.root `ls ${root_file_base}_???.root`



