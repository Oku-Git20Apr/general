#!/bin/bash

function list_recursive()
{
	#defined as a local variable
	local filepath=$1

	echo "$filepath"

	#in the case of directory
	if [-d "$filepath"]; then
		local fname
		for fname in $(ls "$filepath")
		do
			#echo "${filepath}/${fname}"
			list_recursive "${filepath}/${fname}"
		done
	fi
