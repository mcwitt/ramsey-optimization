#!/bin/sh

R=$1
S=$2
NV=$3
seed=$4

cmd="./pt_$R-$S-${NV}.out T.txt 500000 $seed"
prefix=$R-$S-${NV}_$seed
dest=zeros

$cmd > $prefix.tmp

if [ $? -eq 0 ]
then
    mv $prefix.tmp $dest/$prefix.log
    mv $prefix.graph $dest
else
    rm $prefix.tmp
    rm $prefix.graph
fi


