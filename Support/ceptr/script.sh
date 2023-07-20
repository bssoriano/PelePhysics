#!/bin/bash
#SBATCH --account=hpacf
#SBATCH --time=00:59:00
#SBATCH --partition=debug
#SBATCH --job-name=lean_noj
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --qos=high

module purge
module load conda
conda activate /projects/hpacf/mhassana/QSSPaper/conda_envs/mech_format
source /projects/hpacf/mhassana/QSSPaper/PelePhysics/setup.sh

mech_name=dodecane_lu_qss

rm -rf ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/$mech_name_*
cd ${PELE_PHYSICS_HOME}/Support/ceptr
rm -rf log_expr
${poetry_bin} run python generate_multi_qss.py -mf ${PELE_PHYSICS_HOME}/Support/Mechanism/Models -tm $mech_name -nm 100 


for dir in ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/$mech_name_*/; do
    srun ${poetry_bin} run qssa -f ${dir}/skeletal.yaml -n ${dir}/non_qssa_list.yaml
    srun ${poetry_bin} run convert -f ${dir}/qssa.yaml --qss_format_input ${dir}/qssa_input.toml --qss_symbolic_jacobian
done

