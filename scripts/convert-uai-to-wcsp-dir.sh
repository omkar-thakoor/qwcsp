#!/bin/bash

if [[ -d $1 ]]; then

    FILES="$1/*.uai"
    echo "${a}"
    for f in $FILES
    do
	filepath_no_ext="${f%.*}"
	filepath="${filepath_no_ext}.wcsp"
        echo "Processing $f"
	python3 convert-uai-to-wcsp.py < $f > "${filepath}"
    done

else
    echo "$1 is not a valid directory"
    exit 1
fi

