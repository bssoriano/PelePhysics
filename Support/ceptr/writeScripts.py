import os
import argparse
import random

parser = argparse.ArgumentParser(description="Generate multiple mechs")
parser.add_argument(
    "-mf",
    "--mech_folder",
    type=str,
    metavar="",
    required=False,
    help="Folder that contains the mechanisms",
    default=None,
)
parser.add_argument(
    "-tm",
    "--template_mech",
    type=str,
    metavar="",
    required=False,
    help="Template mechanism to copy",
    default="dodecane_lu_qss",
)
parser.add_argument(
    "-nm",
    "--number_mechs",
    type=int,
    metavar="",
    required=False,
    help="Number of mechanisms to generate",
    default=2,
)
args = parser.parse_args()


f = open("script.sh", "w+")


f.write("#!/bin/bash\n")
f.write("#SBATCH --account=hpacf\n")
f.write("#SBATCH --time=07:59:00\n")
f.write("##SBATCH --partition=debug\n")
f.write("#SBATCH --job-name=lean_noj\n")
f.write("#SBATCH --mail-type=BEGIN,END,FAIL\n")
f.write("#SBATCH --mail-user=malik.hassanaly@nrel.gov\n")
f.write("#SBATCH --nodes=1\n")
f.write("#SBATCH --ntasks-per-node=1\n")
f.write("##SBATCH --qos=high\n\n")

f.write("module purge\n")
f.write("module load conda\n")
f.write("conda activate /projects/hpacf/mhassana/QSSPaper/conda_envs/mech_format\n")
f.write("source /projects/hpacf/mhassana/QSSPaper/PelePhysics/setup.sh\n")

f.write("cd ${PELE_PHYSICS_HOME}/Support/ceptr\n\n")

for im in range(args.number_mechs):
    f.write("srun ${poetry_bin} run qssa -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/"+args.template_mech+f"_{im}"+"/skeletal.yaml -n ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/"+args.template_mech+f"_{im}"+"/non_qssa_list.yaml\n")
    
    f.write("srun ${poetry_bin} run convert -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/"+args.template_mech+f"_{im}"+"/qssa.yaml --qss_format_input ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/"+args.template_mech+f"_{im}"+"/qssa_input.toml --qss_symbolic_jacobian\n\n")

