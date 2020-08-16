#!/bin/bash

usage()
{
	#get filename of the shell script
	local script_name=$(basename "$0")

	#help
	cat << END
usage: $script_name PATTERN [PATH] [NAME_PATTERN]
find file in current directory recursively, and print lines which match PATTERN.

	PATH		find file in PATH directory, instead of current directory
	NAME_PATTERN	specify name pattern to find file

Examples:
	$script_name return
	$script_name return ~ '*.txt'
END
}

#if no argument
if [ "$#" -eq 0 ]; then
	usage
	exit 1
fi

pattern=$1
directory=$2
name=$3

#2nd parameter is a starting directory
if [ -z "$directory" ]; then
	directory='.' #current dir
fi

#3rd parameter is a searched filetype 
if [ -z "$name" ]; then
	name='*' #any type
fi

#if no directory
if [ ! -d "$directory" ]; then
	echo "$0: ${directory}: No such directory" 1>&2
	exit 2
fi

find "$directory" -type f -name "$name" | xargs grep -nH "$pattern"

