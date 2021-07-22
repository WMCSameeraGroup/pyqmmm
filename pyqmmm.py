#!/usr/bin/env python3

import sys, os, copy
import utils.FileIO as fio
from utils.Periodictable import periodic_table

print("pyQMMM: Running...")

# 1. prepare the environment
# get system variables
script_path = os.path.dirname(sys.argv[0])
cwd = os.getcwd()

# read input.key file
param_path, fix_ein_path, g16_scratch, tinker_path = fio.read_keyfile("input.key")

# find *.EIn file (default g16_scratch=./)
gaussian_ein = fio.getEin_file(ein_location=g16_scratch)

# run if theres a bash file available
if fix_ein_path:
    fio.sh(f"sh {fix_ein_path} {gaussian_ein}")


# 2. read gaussian .EIn file
# print(">Reading Gaussian .EIn file ...")
data = fio.read_file(f"{gaussian_ein}")
n_atoms, coordinates, connectivity, eOu_setting = fio.clean_ein(data)

# remove partial connectivity from tinker connectivity data
connectivity_cleaned = copy.deepcopy(connectivity)
fio.clean_connectivity(connectivity_cleaned)
# fio.remove_files([f"{gaussian_ein}"])


# 3. generate tinker file
# print(">Generating Tinker file ...")
atom_types = fio.read_atom_types(f"{cwd}/atomtypes.dat")

fio.write_txyz(n_atoms, coordinates, connectivity_cleaned, 
               periodic_table=periodic_table, atom_types=atom_types)


# 4. run tinker
# print(">Running Tinker ...")
# forcefield = fio.extract_forcefield()
fio.run_tinker(eOu_setting=eOu_setting, tinker_path=tinker_path)


# 5. construct gaussian Eou file
gradients, dderivatives, ghes = None, None, None
# print(">Generating Gaussian Eou ...")

if eOu_setting >= 0:
    energy, dipole = fio.process_tinker_epout("input")
if eOu_setting >= 1:
    gradients = fio.process_tinker_gout("input", n_atoms)
if eOu_setting == 2:
    # polarizability -> hard coded to 0
    # dipole derivatives -> hard coded to 0
    dderivatives = fio.extract_dipole_derivatives(n_atoms)    
    
    tinker_hessian = fio.extract_hessian("input.hes")      
    matrix = fio.tinker_hessian_matrix(tinker_hessian)      
    index_list = fio.create_indexes(matrix)
    ghes = fio.rearrange_hessian(index_list, tinker_hessian)
    
fio.write_gauEou(g16_scratch=g16_scratch, eOu_setting=eOu_setting, filename=gaussian_ein[:-4], 
                 natoms=n_atoms, energy=energy, dipole=dipole, 
                 gradients=gradients, dipole_derivatives=dderivatives, 
                 gaussian_hessian=ghes)


# 6. remove temporary files and the EIn file
files_to_remove = ["*.EIn", "*.Ein", "*.ein", "*.epout", "*.gout", "*.hes", "*.hesout", "*.xyz"]

# clear scratch folder and the cwd
fio.remove_files(location=g16_scratch, files=files_to_remove)
fio.remove_files(location=cwd, files=files_to_remove)

print("pyQMMM: Data recovery completed!")
