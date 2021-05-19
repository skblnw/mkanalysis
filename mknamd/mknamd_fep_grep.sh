#!/bin/bash

if [[ "$#" -ne 1 ]]; then
    echo ">Usage: mknamd_fep_grep alchemy.fepout"
    exit 0
fi

grep "^#Free energy change" $1 | awk '{print NR" "$12" "$19}' > fepout
