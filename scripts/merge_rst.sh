#!/bin/bash
ID=$1
istr=$2
iend=$3
EXE=../pyathena/merge_rst.py
for (( i=$istr; i<=$iend; i++ ))
do
    num=$(printf "%04d" "$i")
    echo python $EXE -f /u/ckim14/$ID/rst/$ID.$num.rst -d /u/ckim14/$ID/ -i $ID
done
