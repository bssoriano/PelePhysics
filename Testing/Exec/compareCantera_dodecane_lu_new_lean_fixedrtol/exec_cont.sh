#!/bin/bash

rm Pele3d.llvm.ex
rm PPreaction.txt 
rm timing.txt 
rm err.txt

#rm -r aj*
#rm -r gmres*
#rm -r gmresTyp*
#rm -r denseDir*
 
# Find first eyypecutable name available
#make COMP=llvm TPL
make COMP=llvm -j 8
execname=`find . -name "Pele*.ex" | head -1`

#if [ -z "$execname" ]
#then
#    make COMP=llvm -j 8
#    execname=`find . -name "Pele*.ex" | head -1`
#fi


if [[ -f "PPreaction.txt" ]]; then
    rm PPreaction.txt
fi
if [[ -f "log" ]]; then
    rm log
fi
#$execname inputs/inputs_ref

#declare -a ndtArray=(100 500 1000 5000 10000 50000 100000 400000 800000 3200000 10000000 25600000)
declare -a ndtArray=(10000000)
declare -a tolArray=(1e-12 1e-14)
#declare -a ndtArray=(1000 10000 400000)
#declare -a tolArray=(1e-6 1e-8 1e-10)

for ndt in ${ndtArray[@]}; do
    for tol in ${tolArray[@]}; do
        foldAJ_otpt=aj_${tol}_${ndt}
        foldGMRES_otpt=gmres_${tol}_${ndt}
        foldGMRESTyp_otpt=gmresTyp_${tol}_${ndt}
        foldDenseDir_otpt=denseDir_${tol}_${ndt}
        for irep in {1..3}
        do
            $execname inputs/inputsAJ ode.ndt=$ndt ode.atol=$tol hr.outputFolderHR=./${foldAJ_otpt}
            $execname inputs/inputs_gmres ode.ndt=$ndt ode.atol=$tol hr.outputFolderHR=./${foldGMRES_otpt}
            $execname inputs/inputs_gmres_typ ode.ndt=$ndt ode.atol=$tol hr.outputFolderHR=./${foldGMRESTyp_otpt}
            $execname inputs/inputs_dense_direct ode.ndt=$ndt ode.atol=$tol hr.outputFolderHR=./${foldDenseDir_otpt}
        done
        python computeError.py -f ${foldAJ_otpt}/PPreaction.txt -l ${foldAJ_otpt}/err.txt
        python computeError.py -f ${foldGMRES_otpt}/PPreaction.txt -l ${foldGMRES_otpt}/err.txt
        python computeError.py -f ${foldGMRESTyp_otpt}/PPreaction.txt -l ${foldGMRESTyp_otpt}/err.txt
        python computeError.py -f ${foldDenseDir_otpt}/PPreaction.txt -l ${foldDenseDir_otpt}/err.txt
    done
done


#$execname inputs/inputs_dense_direct
#$execname inputs/inputs_gmres
#$execname inputs/inputs_gmres_typ ode.ndt=800000
#$execname inputs/inputs_gmres_typ2



#cd compareCantera/canteraSim
#python homoReact_dodecane_lu.py
#cd ../..
##conda deactivate
#
#cd compareCantera
#python compareResults.py
