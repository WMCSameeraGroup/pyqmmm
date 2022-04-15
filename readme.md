PyQMMM
===

<!-- <img src="figs/PyQMMM.png" alt="drawing" width="400"/> -->

Python interface for subtractive QM/MM calculations with AMOEBA polarizable force field (PyQM/MM) interfaces between Gaussian16 and Tinker programs. PyQM/MM does not require third-party Python packages and therefore works with the default Python installation in the system (version 3.6 or higher).  
Once Gaussian executes an ONIOM(QM:MM) calculation with the keyword `external`,  PyQM/MM converts the Gaussian16 input data (`*.EIn`) into Tinker compatible inputs (`*.xyz`) and then executes the Tinker calculation. Then, the MM data (`*.epout`, `*.gout`, `*.hes`) is recovered to construct a Gaussian16 compatible data source file (`*.EOu`) to continue the ONIOM(QM:MM) calculation. The incorporation of Tinker tools enables users to use a wide range of MM force fields such as AMOEBA, MMFF, MM3, OPLS, CHARMM, AMBER, etc. 

## Installation

Download the program code from <https://github.com/WMCSameeraGroup/pyqmmm> and extract it. (we assume that you have extracted to " `~/`")

```bash
cd ~/pyqmmm
# set permissions
chmod +x pyqmmm.py
```

## Example Calculation

PyQM/MM requires three input files and the parameter file for Tinker. 

1. `*.com`: Gaussian input file.
2. `input.key`: Tinker keywords file. 
3. `atomtypes.dat`: stores Gaussian and Tinker atom types.

### Gaussian16 input
Gaussian uses the `external` keyword to execute external programs for ONIOM(QM:MM) calculations. 

```bash
#p opt(cartesian,maxcyc=100,nomicro) freq=noraman nosymm oniom(wb97xd/6-31G*:external="/home/user/pyqmmmm/pyqmmm.py") geom=connectivity
```

### `input.key` 
Tinker reads  `input.key` to locate the MM parameter file of Tinker (e.g. amoeba09.prm). Also, `g16_scratch` and `tinker_path` would be defined. 

```bash
parameters /path/to/parameters/amoeba09.prm

# Define `g16_scratch` and `tinker_path`.
#----------------------------------------
g16_scratch /path/to/scratch/scr-water
tinker_path /home/useer/apps/tinker
```

### `atomtypes.dat` 

Tinker program needs MM atom types to calculate MM potential energy and derivatives. Therefore, Gaussian16 atom types must be converted to the corresponding Tinker atom types. The `atomtypes.dat` file contains both Gaussian16 and Tinker atom types.

```bash
# file: atomtypes.dat
# atomicNo   element   g16AtomType   tinkerAtomType   atomDescription
  8           O         OW            36                "Water O"
  1           H         HW            37                "Water H"
```

For instance, oxygen in water molecule is identified as `OW` in Gaussian16 and atomtype `36` in amoeba09 force field.

```bash
# file: waterbox.com
...
 O-OW--0.834000  -1   14.23459000   -1.56102400    2.67321600 L
 H-HW-0.417000   -1   13.70916000   -1.77735100    3.45832800 L
 H-HW-0.417000   -1   14.22399200   -0.58589000    2.60858200 L
...
```

```bash
# file: amoeba09.prm
...
atom         36   34    O     "Water O"            8    15.999    2
atom         37   35    H     "Water H"            1     1.008    1
...
```

Then, Gaussian16 calculation can be submitted. 

Submission script:

```bash
# file: run.csh
#!/bin/csh

set INP=waterbox.com
set SCR=wat
setenv GAUSS_SCRDIR /path/to/scratch/folder/${SCR}

mkdir -p $GAUSS_SCRDIR
g16 < ${INP} > ${INP:r}.log
rm -rf $GAUSS_SCRDIR
```

See `./examples/water` for complete set of files.
