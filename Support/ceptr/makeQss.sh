bash disableAllProfile.sh

poetry run qssa -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/skeletal.yaml -n ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/non_qssa_list.yaml

poetry run convert -f ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/dodecane_lu_qss/qssa.yaml --hformat gpu --remove_1 --remove_pow --remove_pow10  --min_op_count 3 --recursive_op_count --store_in_jacobian --round_decimals --recycle_cse
