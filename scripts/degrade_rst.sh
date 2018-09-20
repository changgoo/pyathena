#!/bin/bash
ID=$1
ID2=$2
i=$3
EXE=../pyathena/degrade_rst.py
num=$(printf "%04d" "$i")
echo python $EXE -f /u/ckim14/$ID/$ID.$num.rst -d /u/ckim14/$ID/rst/ -i $ID2
