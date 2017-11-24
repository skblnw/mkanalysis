#!/bin/bash

# 1. Check for names of DCDs
DCDNAME=`grep "mol addfile" vm_measure-findall.tcl | grep "dcd waitfor" | awk '{print $3}' | awk -F '/' '{print $3}' | awk -F '.' '{print $1}'`
echo -e "DCDNAME for vm_measure-findall.tcl is [ $DCDNAME ]"
DCDNAME=`grep "mol addfile" template-vm-measure-pair-tcl | grep "dcd waitfor" | awk '{print $3}' | awk -F '/' '{print $3}' | awk -F '.' '{print $1}'`
echo -e "DCDNAME for template-vm-measure-pair-tcl is [ $DCDNAME ]"

# 2. Check for exist of following analysis programs
prgrams="make_vmdcalc_pair"
for prog in $programs; do
    if ! type $prog; then
        echo "Command '$prog' does not exist."
        exit 1
    fi
done