#!/bin/bash

rm Pele3d.llvm.ex
rm PPreaction.txt 
rm timing.txt 
rm err.txt
rm -r aj*
rm -r gmres*
rm -r gmresTyp*
rm -r denseDir*
 
# Find first eyypecutable name available
#make COMP=llvm TPL
make COMP=llvm -j 8
execname=`find . -name "Pele*.ex" | head -1`

if [ -z "$execname" ]
then
    make COMP=llvm -j 8
    execname=`find . -name "Pele*.ex" | head -1`
fi


if [[ -f "PPreaction.txt" ]]; then
    rm PPreaction.txt
fi
if [[ -f "log" ]]; then
    rm log
fi
$execname inputs/inputs_ref

