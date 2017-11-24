#!/bin/bash

echo "mknamd> Prepare you analysis scripts for NAMD"

SRC_DIR=$HOME/scripts/mkanalysis/Analysis

mkdir -p make_analysis
cd make_analysis
rsync -a $SRC_DIR/combine_dcd .
