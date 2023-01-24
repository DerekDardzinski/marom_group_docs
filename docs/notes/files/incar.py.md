# incar.py
The `incar.py` can be used to generate the INCAR for all of your calculations. To use it, copy the code block below and put it in `~/bin/incar.py` then run:

```bash
chmod +x ~/bin/incar.py
```

Lastly you will need to place your python interpreter path at the top of the file in order for it to work properly as an executible. You can find your path by running the command `which python`:

```bash
(base) usrname(cori) usrname/bin $ which python
/global/homes/d/ddardzin/.local/miniconda3/bin/python
```

Additionally, you will need to add the path to your `potpaw_PBE` file which holds the VASP pseudopotential information.

```python
#!/<PATH TO YOUR PYTHON INTERPETER>
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar
import argparse
import numpy as np
import os

POTCAR_PATH = '<PATH TO potpaw_PBE FOLDER>'

incar_file = """
# general (start)
ALGO = Fast     # Mixture of Davidson and RMM-DIIS algos
PREC = N        # Normal precision
EDIFF = ediff    # Convergence criteria for electronic converge
NELM = 500      # Max number of electronic steps
ENCUT = encut     # Cut off energy
LASPH = True    # Include non-spherical contributions from gradient corrections
NBANDS = nbands    # Number of bands to include in the calculation
BMIX = 3        # Mixing parameter for convergence
AMIN = 0.01     # Mixing parameter for convergence 
SIGMA = 0.05    # Width of smearing in eV

# parallelization
KPAR = kpar        # The number of k-points to be treated in parallel
NCORE = ncore        # The number of bands to be treated in parallel
# general (end)

# scf (start)
ICHARG = 2      # Generate CHG* from a superposition of atomic charge densities
ISMEAR = 0      # Fermi smearing
LCHARG = lcharg   # Write the CHG* files
LWAVE = False   # Does not write the WAVECAR
LREAL = Auto    # Automatically chooses real/reciprocal space for projections
# scf (end)

# opt (start)
ICHARG = 2      # Generate CHG* from a superposition of atomic charge densities
ISMEAR = 0      # Fermi smearing
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR
IBRION = 2      # Ionic relaxation
NSW = 50        # Maximum of 50 ionic steps
# opt (end)

# band (start)
ICHARG = 11     # Calculate eigenvalues from preconverged CHGCAR
ISMEAR = 0      # Fermi smearing
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR files (True for unfolding)
LORBIT = 11     # Projected data (lm-decomposed PROCAR)
# band (end)

# dos (start)
ICHARG = 11     # Calculate eigenvalues from preconverged CHGCAR
ISMEAR = -5     # Tetrahedron method with Blochl corrections
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR files 
LORBIT = 11     # Projected data (lm-decomposed PROCAR)
NEDOS = 3001    # 3001 points are sampled for the DOS
EMIN = emin     # Minimum energy for the DOS plot
EMAX = emax     # Maximum energy for the DOS plot
# dos (end)

# stm (start)
ICHARG = 11     # Calculate eigenvalues from preconverged CHGCAR
ISMEAR = -5     # Tetrahedron method with Blochl corrections
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR file
LPARD = True    # Writes the PARCHG files
NBMOD = -3      # Makes the energy be with respect to the Fermi-energy
EINT = 0 1.5    # Energy range for calculating the PARCHG (bias)
# stm (end)

# wave (start)
ICHARG = 11     # Calculate eigenvalues from preconverged CHGCAR
ISMEAR = -5     # Tetrahedron method with Blochl corrections
LCHARG = False  # Does not write the CHG* files
LWAVE = False   # Does not write the WAVECAR file
LPARD = True    # Writes the PARCHG files
KPUSE = 1       # Determines which k-point(s) to use
IBAND = 255     # Determines which band(s) to use
# wave (end)

# dftu (start)
LDAU = True     # Determines if DFT+U is used
LDAUTYPE = 2    # Dudarev formulation
LDAUL = ldaul   # l-quantum number to apply the U-value on (-1 turns it off)
LDAUU = ldauu   # Effective U-value for each species
LDAUJ = ldauj   # J-value (Always zero for Dudarev method)
LMAXMIX = 4     # Max l-quantum number for charge density mixing
# dftu (end)

# soc (start)
LSORBIT = True  # Turn on spin-orbit coupling
MAGMOM = magmom # Set the magnetic moment for each atom (3 for each atom)
# soc (end)

# sp (start)
MAGMOM = magmom # Set the magnetic moment for each atom (1 for each atom) 
ISPIN = 2       # Spin-polarized calculation
# sp (end)

# slab (start)
IDIPOL = 3      # Calculates the dipole along the z-axis
LREAL = Auto    # Automatically chooses rea/reciprocal space for projections
DIPOL = dipol   # Defining the location of the center of the dipole moment (center of mass)
LDIPOL = True   # Adding dipole corrections
# slab (end)

# hse (start)
LHFCALC = True  # Determines if a hybrid functional is used
HFSCREEN = 0.2  # Range-separation parameter
AEXX = 0.25     # Fraction of exact exchange to be used
PRECFOCK = Fast # Increases the speed of HSE Calculations
# hse (end)

# vdw (start)
IVDW = 20           # Tkatchenko-Scheffler Method
LVDW_EWALD = True   # Use Ewald's Summation
# vdw (end)

# efield (start)
EFIELD = efield
# efield (end)
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description='param')
    parser.add_argument(
        '-s',
        '--scf',
        help="Generates an INCAR for a SCF calculation",
        action='store_true'
    )
    parser.add_argument(
        '-b',
        '--band',
        help="Generates an INCAR for a band structure calculation",
        action='store_true'
    )
    parser.add_argument(
        '-d',
        '--dos',
        help="Generates an INCAR for a density of states calculation",
        action='store_true'
    )
    parser.add_argument(
        '-t',
        '--stm',
        help="Generates an INCAR for a STM calculation",
        action='store_true'
    )
    parser.add_argument(
        '-w',
        '--wave',
        help="Generates an INCAR for a wave function calculation",
        action='store_true'
    )
    parser.add_argument(
        '-o',
        '--opt',
        help="Generates an INCAR for an opt calculation",
        action='store_true'
    )
    parser.add_argument(
        '-u',
        '--dftu',
        help="Adds options for DFT+U into the INCAR",
        action='store_true'
    )
    parser.add_argument(
        '-c',
        '--soc',
        help="Adds options for spin-orbit coupling into the INCAR",
        action='store_true'
    )
    parser.add_argument(
        '-p',
        '--sp',
        help="Adds options for a spin-polarized calculation into the INCAR",
        action='store_true'
    )
    parser.add_argument(
        '-l',
        '--slab',
        help="Adds options for a slab calculation into the INCAR",
        action='store_true'
    )
    parser.add_argument(
        '-e',
        '--hse',
        help="Adds options for an HSE calculation into the INCAR",
        action='store_true'
    )
    parser.add_argument(
        '-v',
        '--vdw',
        help="Adds options for including Van der Waals interactions into the INCAR",
        action='store_true'
    )
    parser.add_argument(
        '-f',
        '--efield',
        help="Adds options for including an electric field into the INCAR",
        action='store_true'
    )
    parser.add_argument(
        '--output',
        help="Output file name",
        default='INCAR',
        type=str,
    )
    parser.add_argument(
        '--poscar',
        help="POSCAR used for setting the MAGMOM tag for SOC Calculations",
        default='POSCAR',
        type=str,
    )
    parser.add_argument(
        '-m',
        '--magmom',
        help="Magnetic moments for setting the MAGMOM tag given the POSCAR. The MAGMOMs should be in the same order as shown in the POSCAR (e.g. In Sb Fe --> --magmom 0 0 2.18)",
        type=float,
        nargs='+',
    )
    parser.add_argument(
        '--encut',
        help="Used to set the value of the ENCUT",
        default=400,
        type=int,
    )
    parser.add_argument(
        '--kpar',
        help="Used to set the value of KPAR",
        default=8,
        type=int,
    )
    parser.add_argument(
        '--ncore',
        help="Used to set the value of NCORE",
        default=8,
        type=int,
    )
    parser.add_argument(
        '--ldaul',
        help="L values for setting the LDAUL tag given the POSCAR. The Ls should be in the same order as shown in the POSCAR (e.g. In As --> --ldaul 1 1)",
        type=int,
        nargs='+',
    )
    parser.add_argument(
        '--ldauu',
        help="U values for setting the LDAUU tag given the POSCAR. The Ls should be in the same order as shown in the POSCAR (e.g. In As --> --ldauu -0.5 -7.5)",
        type=float,
        nargs='+',
    )
    parser.add_argument(
        '--lcharg',
        help="Option to turn of LCHARG",
        default=1,
        type=int,
    )
    parser.add_argument(
        '--ediff',
        help="Option to set EDIFF",
        default='1e-5',
        type=str,
    )
    parser.add_argument(
        '--addbands',
        help="Used to add an additonal amount of bands to the calculation",
        default=0,
        type=int,
    )
    return parser.parse_args()


args = parse_arguments()
output = args.output
poscar = args.poscar
magmom = args.magmom
ldaul = args.ldaul
ldauu = args.ldauu
lcharg = args.lcharg
ediff = args.ediff
encut = args.encut
kpar = args.kpar
ncore = args.ncore
addbands = args.addbands

keys_array = list(args.__dict__.keys())[:-11]
values_array = list(args.__dict__.values())[:-11]
selected_options = [k for k, v in zip(keys_array, values_array) if v]

if args.hse and args.scf:
    selected_options.remove('hse')

if args.stm and args.scf:
    selected_options.remove('stm')

if args.wave and args.scf:
    selected_options.remove('wave')

incar = incar_file.split('\n')[1:-1]

start_inds = []
end_inds = []
keys = []

for i, line in enumerate(incar):
    if '#' and 'start' in line:
        start_inds.append(i)
        keys.append(line.split()[1])
    if '#' and 'end' in line:
        end_inds.append(i)

incar_dict = {
    k: incar[s:e] for s, e, k in zip(start_inds, end_inds, keys)
}

selected_incar = incar_dict['general']

for k in selected_options:
    selected_incar.extend([''] + incar_dict[k])

selected_incar.extend([''])

if args.hse and not args.scf:
    selected_incar[1] = 'ALGO = D        # Damped velocity friction algorithm'

new_incar = '\n'.join(selected_incar)
new_incar = new_incar.replace('(start)', '')

if args.hse and not args.scf:
    new_incar = new_incar.replace('ICHARG = 11', 'ICHARG = 0 ')

if args.hse and args.scf:
    new_incar = new_incar.replace('LWAVE = False', 'LWAVE = True ')

if args.stm or args.wave:
    if not args.scf:
        new_incar = new_incar.replace('KPAR = 4', 'KPAR = 1')
        new_incar = new_incar.replace('NPAR = 7', 'NPAR = 1')
    if args.scf:
        new_incar = new_incar.replace('LWAVE = False', 'LWAVE = True ')

new_incar = new_incar.replace('ENCUT = encut', f'ENCUT = {encut}')
new_incar = new_incar.replace('KPAR = kpar', f'KPAR = {kpar}')
new_incar = new_incar.replace('NCORE = ncore', f'NCORE = {ncore}')
new_incar = new_incar.replace('LCHARG = lcharg', f'LCHARG = {str(bool(lcharg))}')
new_incar = new_incar.replace('EDIFF = ediff', f'EDIFF = {ediff}')

if args.dftu:
    if args.ldaul is not None:
        new_incar = new_incar.replace('LDAUL = ldaul', f'LDAUL = {" ".join([str(l) for l in ldaul])}')
        new_incar = new_incar.replace('LDAUJ = ldauj', f'LDAUJ = {" ".join([str(0) for _ in ldaul])}')

    if args.ldauu is not None:
        new_incar = new_incar.replace('LDAUU = ldauu', f'LDAUU = {" ".join([str(u) for u in ldauu])}')
        new_incar = new_incar.replace('LDAUJ = ldauj', f'LDAUJ = {" ".join([str(0) for _ in ldauu])}')

if os.path.isfile(poscar):
    p = Poscar.from_file(poscar, check_for_POTCAR=False)
    elements = p.site_symbols
    natoms = p.natoms
    nelect = []

    for element in elements:
        with open(os.path.join(POTCAR_PATH, element, 'POTCAR'), 'r') as potcar:
            potcar.readline()
            zval = float(potcar.readline().split()[0])
            nelect.append(zval)

    natoms = np.array(natoms)
    nelect = np.array(nelect)

    nbands_1 = int((natoms * nelect).sum() + natoms.sum()) + addbands
    nbands_2 = int(1.2 * (natoms * nelect).sum()) + addbands

    if args.soc:
        nbands = int(np.max([nbands_1, nbands_2]))
    else:
        nbands = int(0.5*np.max([nbands_1, nbands_2]))

    new_incar = new_incar.replace('NBANDS = nbands', f'NBANDS = {nbands}')

    if args.slab:
        weights = [s.species.weight for s in p.structure]
        center_of_mass = np.average(p.structure.frac_coords, weights=weights, axis=0)

        com = []
        for i in center_of_mass:
            com.append(f'{round(i,5)}')

        new_incar = new_incar.replace('DIPOL = dipol', f'DIPOL = {" ".join(com)}')

    if args.soc:
        if args.magmom is None:
            new_incar = new_incar.replace('MAGMOM = magmom', f'MAGMOM = {int(3 * natoms.sum())}*0')
        else:
            symbols = p.site_symbols
            natoms = p.natoms

            mags = []
            for i, s in enumerate(symbols):
                temp_mags = []
                for _ in range(natoms[i]):
                    temp_mags.append(f'0 0 {args.magmom[i]}')

                mags.extend(temp_mags)

            new_incar = new_incar.replace('MAGMOM = magmom', f'MAGMOM = {" ".join(mags)}')

    if args.sp:
        if args.magmom is None:
            new_incar = new_incar.replace('MAGMOM = magmom', f'MAGMOM = {int(natoms.sum())}*0')
        else:
            symbols = p.site_symbols
            natoms = p.natoms

            mags = []
            for i, s in enumerate(symbols):
                mags.append(f'{int(natoms[i])}*{args.magmom[i]}')

            new_incar = new_incar.replace('MAGMOM = magmom', f'MAGMOM = {" ".join(mags)}')


else:
    print('The POSCAR file does not exist, please make sure the file name is correct.')


with open(output, 'w') as out:
    out.write(new_incar)
```

