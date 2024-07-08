#!/usr/bin/env python3

import sys


def caldb_hp_masked_or_not(file_path, target_string):
    with open(file_path, "r") as file:
        lines = file.readlines()
        if len(lines) >= 9 and target_string in lines[8]:
            return True
        else:
            return False


L1_dataset = str(sys.argv[1])

drivermodule_param_file = "UVIT_DriverModule.par"
target_string = "frm05274"

mask_needed = int(L1_dataset[35:40]) >= 5274
hp_masked = caldb_hp_masked_or_not(drivermodule_param_file, target_string)

if hp_masked == mask_needed:
    print("CALDB configuration is correct.")
else:
    raise RuntimeError("STOP! Check CALDB configuration.")
