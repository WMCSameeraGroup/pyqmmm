import sys, os, shlex, subprocess, glob
import re 
# import itertools

SCRIPT_PATH = os.path.dirname(__file__)
TINKER_TIMEOUT = float(500)  # ? increase this later?

# unit conversion
BOHR2ANG = float(0.529177)  # Coordinates: gaussian16 (Bohr) -> tinker (Angstrom)
KCALPERMOL2HARTREE = float(0.001593601)  # Energy: tinker (kcal mol-1) -> gaussian16 (Hartree)
DEBYE2BOHRELEC = float(0.393430343)  # Dipole moment: tinker (Debye) -> gaussian16 (Bohr-electron)
KCALPERMOLANG2HARTREEBHOR  = float(0.000843297)  # Gradient: tinker (kcal mol-1/Å) -> gaussian (Hartree/Bohr)
KCALPERMOLANG22HARTREEBHOR2 = float(0.000446254)  # Hessian: tinker (kcal mol-1/Å^2) -> gaussian (Hartree/Bohr^2)


def abnormal_termination():
    print(">Terminated with an error!")
    sys.exit()


def read_file(file_name):
    """(str -> list)
    Read gaussian .Ein file to a list.
    Does not strip() or split()

    Args:
        file_name (str): path to .dat file
    """
    try:
        with open(file_name, "r") as f:
            file_data = f.read().splitlines()
    except:
        print(f"ERROR: Can not open {file_name} file!")
        print("(file doesn't exist or may be in a different path)")
        abnormal_termination()

    return(file_data)


def clean_ein(data):
    """(list -> nested list)
    Clean gaussian .Ein file and split it to
     - number of atoms: nAtoms
     - xyz coordinates: coordinates
     - atom connectivity: connectivity

    Args:
        data (list): gaussian file
    """
    line0 = data[0].split() # read first line and split
    nAtoms = int(line0[0])  # number of atoms -> line [0]
    eOu_setting = int(line0[1]) # what to write in EOu file -> line [1]
    coordinates = []  # coordinate data -> line [1:nAtoms+1]
    # connectivity data -> line [(nAtoms+1):(nAtoms+1)+(nAtoms+1)]
    connectivity = []

    # extract coordinates. Line e.g.
    try:
        for line in data[1:(nAtoms+1)]:
            coordinates_line_split = line.split()
            if len(coordinates_line_split) == 6:
                coordinates.append(coordinates_line_split)
            else:
                raise ValueError
        
    except ValueError as error:
        print(">Error in Gaussian .EIn file")
        print(f"in line: {line}")
        abnormal_termination()        

    # extract connectivity
    for line in data[(nAtoms+1):(nAtoms+1)+(nAtoms+1)]:
        connectivity.append(line.split())

    return(nAtoms, coordinates, connectivity, eOu_setting)


def read_atom_types(file_name):
    """(str -> nested list)
    Create atom types dic from a .dat file
     - key: Gaussian MM level atom type
     
    Helps to convert Gaussian MM level atom type(s) to Tinker atom type(s).
    NOTE: atom type depends on the force field.

    Args:
        file_name (str): path to .dat file 
    """
    atom_types = {}
    file_data = read_file(file_name)

    for line in file_data:
        lineSplit = line.strip().split()
        # skip empty lines
        if lineSplit:
            # skip comment lines starts with # (including "# " and "#ABC")
            if not lineSplit[0][0] == "#":
                try:
                    atomicNo = lineSplit[0]
                    element = lineSplit[1]
                    g16AtomType = lineSplit[2]
                    tinkerAtomType = lineSplit[3]
                except IndexError as error:
                    print(">ERROR: Check atomtypes.dat file")
                    abnormal_termination()

                # * this updates/replaces if same key is added
                # * only later value remains in the dict
                atom_types.update({g16AtomType: {
                    "element": element,
                    "atomicNo": int(atomicNo),
                    "tinkerAtomType": int(tinkerAtomType)
                }})
    return(atom_types)


def write_txyz(n_atoms, coordinates, connectivity, periodic_table, atom_types):
    """(int, list... -> NONE)
    Dump Tinker .xyz file

    Args:
        n_atoms (int): total number of atoms in the system
        coordinates (list): xyz coordinates
        connectivity (list): atom connectivity information
        periodic_table (dic): atomic numbers and elements
        atom_types (dic): force field / gaussian atom types  
    """
    # ! writes input.xyz in the working directory
    # TODO: dynamically update force field in the description "(FF param)"
    xyz_description = "pyQMMM generated tinker .xyz (XXX param)"
    with open("input.xyz", "w") as file_out:       
        file_out.write(f"{n_atoms}  {xyz_description}\n")
        for idx in range(n_atoms):
            try:
                column1 = idx+1  # atom index
                column2 = periodic_table[int(coordinates[idx][0])]["elementSymbol"]  # element
                column3 = float(coordinates[idx][1]) * BOHR2ANG  # coordinate x (Angstrom)
                column4 = float(coordinates[idx][2]) * BOHR2ANG  # coordinate y (Angstrom)
                column5 = float(coordinates[idx][3]) * BOHR2ANG  # coordinate z (Angstrom)
                column6 = atom_types[coordinates[idx][5]]["tinkerAtomType"]  # force field atom type
                # atom connectivity
                column7 = " ".join(map(str, connectivity[idx][1::2]))
            
            except KeyError as kerr:
                print(f">ERROR: check .EIn file! can't find atomic number {coordinates[idx][0]} in Periodictable.py")
                print(f"or atomtype {coordinates[idx][5]} in MM force field file")
                abnormal_termination()

            row = " {:d} {:3} {:14.8f} {:14.8f} {:14.8f} {:4}  {} \n".format(
                column1, column2, column3, column4, column5, column6, column7)
            file_out.write(row)


def sh(command, timeout=120):
    """(str -> NONE)
    Execute shell command

    Args:
        command (str): shell command to execute
    """
    shell_arguments = shlex.split(command)
    print(">Running a quick fix...")
    try:
        process = subprocess.Popen(shell_arguments, stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   universal_newlines=True)
        process_output, process_error = process.communicate(timeout)
        process_status = process.wait()
        if process_error:
            print(f">>{process_error}")
            raise OSError
        if process_output:
            print(f">>{process_output}")
    except(OSError, ValueError) as error:
        print(f"ERROR: Can't execute {command}! \nExiting...")
        abnormal_termination()
    except(subprocess.TimeoutExpired) as error:
        print(f"ERROR: Took too long (>{timeout} sec) to process! Terminating...")
        process.terminate()
        abnormal_termination()


def execute_tinker(command, out_file):
    """(str -> NONE)
    Execute shell command

    Args:
        command (str): shell command to execute
    """
    shell_arguments = shlex.split(command)
    try:
        with open(out_file, "w") as fout:
            process = subprocess.Popen(shell_arguments, stdout=fout)
            process_output, process_error = process.communicate(timeout=TINKER_TIMEOUT)
            process_status = process.wait()
    except OSError as error:
        print(">ERROR in Tinker calculation")
        print("can't find executables. Add then in to PATH or specify in the .key file ")
        print("to identify the problem")
        abnormal_termination()
        
    except ValueError as error:
        print(">ERROR in Tinker calculation")
        print("check *.epout  *.gout  *.hes  *.hesout *.xyz files in the working directory")
        print("to identify the problem")
        abnormal_termination()
        
    except(subprocess.TimeoutExpired) as error:
        print(f"ERROR: Tinker took too long (>{TINKER_TIMEOUT} sec) to process! Terminating...")
        print("probable reason: Tinker might not be able to find forcefield parameters. Check")
        print(" - *.key file")
        print(" - atomtypes.dat")
        process.terminate()
        abnormal_termination()


def read_keyfile(file_name):
    """Check if *.key file exsists and read parameters path
    and fix_ein shell script path

    Args:
        file_name (str): full/relative path to *.key file
    """
    param_path   = ""
    fix_ein_path = ""
    tinker_path  = ""
    g16_scratch  = "./"
    
    
    # check if input.key is available if not terminate
    if not os.path.isfile(file_name):
        print(f"ERROR: Can't read {file_name}")
        abnormal_termination()

    keyfile_data = read_file(file_name)
    for line in keyfile_data:
        try:
            line_split = line.strip().split()
            if line_split[0].lower() == "parameters":
                param_path = line_split[1]
            if line_split[0].lower() == "fix_ein":
                fix_ein_path = line_split[1]
            if line_split[0].lower() == "g16_scratch":
                g16_scratch = f"{line_split[1]}/"
            if line_split[0].lower() == "tinker_path":
                tinker_path = f"{line_split[1]}/"
        except:
            pass        
    return(param_path, fix_ein_path, g16_scratch, tinker_path)


def run_tinker(eOu_setting, filename="input", tinker_path=""):
    """Run Tinker executables.
     - analyze
     - testgrad
     - testhess
    Executables must be installed in the target system.
    Force field must be defined in the *.key file
    """         
    if eOu_setting==0:
        # run analyze
        execute_tinker(f"{tinker_path}analyze input.xyz E,M", out_file=f"{filename}.epout")
    elif eOu_setting==1:
        # run analyze and testgrad
        execute_tinker(f"{tinker_path}analyze input.xyz E,M", out_file=f"{filename}.epout")
        execute_tinker(f"{tinker_path}testgrad input.xyz y n 0.1D-04", out_file=f"{filename}.gout")
    elif eOu_setting==2:
        execute_tinker(f"{tinker_path}analyze input.xyz E,M", out_file=f"{filename}.epout")
        execute_tinker(f"{tinker_path}testgrad input.xyz y n 0.1D-04", out_file=f"{filename}.gout")
        execute_tinker(f"{tinker_path}testhess input.xyz y n", out_file=f"{filename}.hesout")
    else:
        print("Someting went wrong. Tinker failed to run.")
        abnormal_termination()


def process_tinker_epout(input_filename="input"):
    """Extract data from Tinker outputs

    Args:
        input_filename (str): Tinker output file name
    """
    energy, dx, dy, dz = None, None, None, None
    
    # read *.epout
    data_tinker_epout = read_file(f"{input_filename}.epout")
    
    # extract energy
    # *line=" Total Potential Energy :   5.9172 Kcal/mole"
    for line in data_tinker_epout:
        if line.startswith(" Total Potential Energy :"):
            # print("Extracting MM Energy:")
            # print(line)
            energy = float(line.split()[4]) * KCALPERMOL2HARTREE
    
    # extract dipole
    # *line=" Dipole X,Y,Z-Components :   0.000   0.000   0.000"
        if line.startswith(" Dipole X,Y,Z-Components :"):
            # print("Extracting MM Dipole:")
            # print(line)
            line_split= line.split()
            dx = float(line_split[3]) * DEBYE2BOHRELEC
            dy = float(line_split[4]) * DEBYE2BOHRELEC
            dz = float(line_split[5]) * DEBYE2BOHRELEC
            break   
    dipole = [dx, dy, dz]
    
    if not (energy and dipole):
        print("Can't extract energy or dipole from Tinker output")
        abnormal_termination()
        
    return(energy, dipole)


def process_tinker_gout(input_filename, n_atoms):
    """Extract data from Tinker outputs

    Args:
        input_filename (str): Tinker file name
        n_atoms (int): number of atoms in the system
    """               
    # read *.gout
    data_tinker_gout = read_file(f"{input_filename}.gout")
    
    # extract gradient
    gradients = []
    """ *
      Type      Atom            > dE/dX     > dE/dY     > dE/dZ      Norm

    Anlyt         1           -17.3537     49.0858     -0.0000       52.0631
    Anlyt         2            -8.9795    -21.3780    -27.0056       35.5943
    ...
    """
    for idx, line in enumerate(data_tinker_gout):
        if line.startswith("  Type      Atom              dE/dX"):
            # print("Extracting MM Gradient:")
            # print(idx, line)
            line_index = idx
            break
    try:
        st =  line_index+2 ; ed = st +  n_atoms
    except:
        print(f"Can't find cartesian gradient breakdown in the {input_filename}.gout file")
        abnormal_termination()
        
    for line in range(st,ed):
        # print(line, data_tinker_gout[line])
        data_gradient_split = data_tinker_gout[line].split()
        dEdx = float(data_gradient_split[2]) * KCALPERMOLANG2HARTREEBHOR
        dEdy = float(data_gradient_split[3]) * KCALPERMOLANG2HARTREEBHOR
        dEdz = float(data_gradient_split[4]) * KCALPERMOLANG2HARTREEBHOR
        gradients.append([dEdx, dEdy, dEdz])
        
    return(gradients)


def getEin_file(ein_location="./"):
    """(empty -> str)
    Automatically detect *.EIn filetype and return the file name. 
    
    Args:
        ein_location (str): location where *.EIn is. Normally Gaussian scratch directory

    Raises:
        Exception: if multiple *.EIn files are present in the working directory.
    """
    gauEin_files = []
    filetypes=["*.Ein", "*.ein", "*.EIn"]
    
    for filetype in filetypes:
        search_pattern = f"{ein_location}/{filetype}"
        gauEin_files.extend(glob.glob(search_pattern))
        
    # filter correct *.Ein file. Especially if there are multiples
    # usual file name is "Gau-1234.Ein"
    # TODO: rectify this later
    # for file in gauEin_files:
    #     if not file[:3] == "Gau":
    #         gauEin_files.remove(file)
            
    try:
        # check if there are multiple *.Ein files, which skipped 
        # through the initial filtration
        if len(gauEin_files) > 1:
            raise Exception()
        gauEin_filename = gauEin_files[0]
    except:
        print("ERROR:", gauEin_files, ein_location) 
        print("   - Can't find *.Ein file or")
        print("   - There are multiple *.Ein files")
        abnormal_termination()
    
    return(gauEin_filename)
    

def remove_files(files, location):
    """(list -> NONE)
    Remove files given as a list.
    
    Args:
        files (list): list of files to be removed
    """
    rm_files_list = []
    for file in files:
        search_pattern = f"{location}/{file}"
        rm_files_list.extend(glob.glob(search_pattern))
        
    for file in rm_files_list:
        os.remove(file) 
    print(f"{len(rm_files_list)} File(s) removed from {location}")



def check_tinker_executables():
    """
    # TODO: Check if Tinker executables are in the $PATH
    """
    pass


# def extract_forcefield():
#     """(str -> str)
#     Extract force field from the key file.
#     """
#     key_file = glob.glob("*.key")
#     if len(key_file) > 1:
#         print(f"WARNING: Multiple *.key files are present. {key_file[0]} selected!")
#     keyfile_data = read_file(key_file[0])

#     try:
#         for line in keyfile_data:
#             if line.startswith("parameters"):
#                 forcefield_path = line.split()[1]
#         return(forcefield_path)
#     except:
#         print("ERROR: Force field path is not defined in the Tinker *.key file")
#         print("Terminating...")
#         abnormal_termination()


def extract_hessian(filename):
    """(text file -> list)
    Extracts hessian from tinker output as it is.    

    Args:
        filename (str): Tinker hessian output file name E.g. inp.hes
    """    
    # print("Extracting Hessian...")
    
    tinker_hessian = []
    hessian_data = read_file(filename) 
    
    iter_hes_data = iter(hessian_data)
    for idx, line in enumerate(iter_hes_data):
        row = []
        if "H" in line:
            # print(idx, line)
            # extract hes data rows
            hes_idx = idx+2
            try:
                while not hessian_data[hes_idx] == "":
                    # print(f"({hes_idx}) {hessian_data[hes_idx]}")
                    row.extend(hessian_data[hes_idx].split())
                    hes_idx += 1
            except IndexError as idx_error:
                # print("Reached to the end...")
                pass
            # print(f"idx={idx}, hes_idx={hes_idx}")
            tinker_hessian.append(row)
            # print(row)
            # TODO: performance enhance with itertools
            # next(itertools.islice(iter_hes_data, hes_idx, None), None) 
    return(tinker_hessian)


def tinker_hessian_matrix(tinker_hessian):
    """nested list -> nested list
    Arrange extracted tinker hessian output to a symmetric matrix 

    Args:
        tinker_hessian (nested list): hessian directly from tinker.
    """
    matrix = []
    diagonal_elements = len(tinker_hessian[0])

    for i in range(1, diagonal_elements):
        temp_list=[]
        for j in range(i+1):
            # print(j)
            temp_list.append(j)
        matrix.append(temp_list)
        # print("")
        
    return(matrix)


def create_indexes(matrix):
    """nested list -> nested list
    Create list index list to convert tinker hessian to gaussian Eou format:
    refer: https://gaussian.com/external/

    Args:
        matrix (nested list): symmetric matrix of tinker hessian
    """
    index_list = []
    for row in matrix:
        reversed_row = row.copy()
        reversed_row[-1] = reversed_row[-1] -1 
        reversed_row.reverse()
        temp_index_list = zip(row, reversed_row)
        index_list.extend(temp_index_list)  
    #add final element (0,n)
    index_list.append((0,len(matrix)))      

    return(index_list)


def rearrange_hessian(index_list, tinker_hessian):
    """nested list, nested list -> nested list
    Generate gaussian Eou hessian from tinker hessian

    Args:
        index_list (nested list): index list to extract and convert tinker hessian
            to gaussian format 
        tinker_hessian (nested list): hessian (symmetric matrix) from tinker
    """      
    # rearrange hessian
    gaussian_hessian = []
    for element in index_list:
        a = element[0]
        b = element[1]
        temp_hes = tinker_hessian[a][b] 
        temp_hes = float(temp_hes) * KCALPERMOLANG22HARTREEBHOR2
        gaussian_hessian.append("{:20.12e}".format(temp_hes))
        # print(temp_hes)
    return(gaussian_hessian)


def extract_polarizability():
    """
    pyQMMM does not extract polarizability information.
    hard-coded to 0
    """
    return(0)


def extract_dipole_derivatives(natoms):
    """
    pyQMMM does not extract dipole_derivatives information.
    hard-coded to 0

    Args:
        natoms (int): number of atoms in the system
    """
    dderivatives = ["{:20.12e}".format(0)] * 9 * natoms
    return(dderivatives)


def write_gauEou(g16_scratch, eOu_setting, filename, natoms, energy, dipole, gradients, gaussian_hessian, dipole_derivatives):
    """(str, int, str, float, list, nested list -> NONE)
    remove old *.EOu files and 
    write new gaussian *.EOu file in the working/scratch directory 

    Args:
        g16_scratch (str):
        eOu_setting (int):
        filename (str): output file name (e.g. g16.eou)
        energy (float): MM total potential energy
        dipole (list): MM dipole
        gradients (nested list): MM gradient
        gaussian_hessian (nested list): hessian
    """
    # remove existing *.Eou files
    remove_files(location=g16_scratch, files=["*.EOu"])
    
    # write new *.EOu file
    with open(f"{filename}.EOu", "w") as fout:
        # write gaussian file
        # print(f">Writing {filename}.EOu")
        
        if eOu_setting >= 0:
            # 1. energy, dipole-moment (xyz)
            fout.write("{:20.12e}{:20.12e}{:20.12e}{:20.12e} \n".format(energy, dipole[0], dipole[1], dipole[2]))
        
        if eOu_setting >= 1:
            # 2. gradient on atom (xyz)
            for line in gradients:
                fout.write("{:20.12e}{:20.12e}{:20.12e} \n".format(line[0], line[1], line[2]))
                  
        if eOu_setting == 2:    
            # 3. polarizability
            # ! hard coded to 0
            pol = 0.0
            fout.write("{:20.12e}{:20.12e}{:20.12e} \n".format(pol, pol, pol))
            fout.write("{:20.12e}{:20.12e}{:20.12e} \n".format(pol, pol, pol))  
                
            # 4. dipole derivatives
            # ! hard coded to 0
            count = 0
            while count < len(dipole_derivatives):
                fout.write(f"{dipole_derivatives[count]}{dipole_derivatives[count+1]}{dipole_derivatives[count+2]} \n")
                count += 3   
                
            # 5. hessian
            count = 0
            while count < len(gaussian_hessian):
                fout.write(f"{gaussian_hessian[count]}{gaussian_hessian[count+1]}{gaussian_hessian[count+2]} \n")
                count += 3
