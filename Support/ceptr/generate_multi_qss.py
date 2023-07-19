import shutil
from distutils.dir_util import copy_tree
from shutil import copy
import toml
import os
import numpy as np
import itertools

import argparse

import random

parser = argparse.ArgumentParser(description="Generate multiple mechs")
parser.add_argument(
    "-mf",
    "--mech_folder",
    type=str,
    metavar="",
    required=True,
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


choices = {}
choices['h_format'] = ['gpu', 'cpu']
choices['remove_1'] = [True, False]
choices['remove_pow'] = [True, False]
choices['remove_pow10'] = [True, False]
choices['round_decimals'] = [True, False]

choices['min_op_count'] = [0]
choices['min_op_count_all'] = list(range(10))
choices['gradual_op_count'] = [True]
choices['remove_single_symbols_cse'] = [True]

choices['store_in_jacobian'] = [True, False]
choices['recycle_cse'] = [True, False]

keys = [key for key in choices]
c = [choices[key] for key in keys]

allComb = list(
    itertools.product(*c)
)
sweep = random.sample(allComb, args.number_mechs)
np.savez("params.npz", val=sweep, name=keys)


template_toml = toml.load(os.path.join(args.mech_folder, args.template_mech, 'qssa_input.toml'))

for im in range(args.number_mechs):
    copy_tree(
        os.path.join(args.mech_folder, args.template_mech),
        os.path.join(args.mech_folder, args.template_mech + f"_{im}")
    )
    new_toml = template_toml.copy()
    for ik, key in enumerate(keys):
        found = False
        for upperKey in new_toml:
            if key in new_toml[upperKey]:
                new_toml[upperKey][key] = sweep[im][ik]
                found = True
            if found:
                break

    file_name = os.path.join(args.mech_folder, args.template_mech + f"_{im}", 'qssa_input.toml')
    with open(file_name, "w") as toml_file:
        toml.dump(new_toml, toml_file)

