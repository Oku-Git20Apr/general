#!/bin/bash

function list_recursive()
{
	#defined as a local variable
	local filepath=$1
	local indent=$2

	#show relative pathname w/ indent
	echo "${indent}${filepath##./}"

	#in the case of directory
	if [ -d "$filepath" ]; then
		local fname
		
		_IFS=$IFS #IFS backup
		#IFS is defined as $'\n'
		IFS='
		'

		for fname in $(ls "$filepath")
		do
			#echo "${filepath}/${fname}"
			list_recursive "${filepath}${fname}" "		$indent"
		done
	
		IFS=$_IFS
	fi

}

list_recursive "$1"
