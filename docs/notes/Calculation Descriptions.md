# Calculation Descriptions
In this section we will discus the inputs you need in order to run all the different types of calculations that we run in the Marom group.

## General Inputs
### INCAR
As stated in [[VASP Basics]], the INCAR contains the set of instructions that tells VASP what type of  calculation to perform. Although we perform many different types of calculations, there are some standard parameters that we always use. Sometimes, we may choose to alter   the following parameters slightly, but for the vast majority of our calculations, these will remain constant.

```txt
ALGO = Fast      # Mixture of Davidson and RMM-DIIS algos
PREC = N         # Normal precision
EDIFF = 1e-5     # Convergence criteria for electronic converge
NELM = 500       # Max number of electronic steps
ENCUT = 350      # Cut off energy
LASPH = True     # Include non-spherical contributions from gradient corrections
NBANDS = 100     # Number of bands to include in the calculation
BMIX = 3         # Mixing parameter for convergence
AMIN = 0.01      # Mixing parameter for convergence
SIGMA = 0.05     # Width of smearing in eV
```

### POSCAR
The POSCAR is the most user-dependent file in VASP and it defines the unit cell and exact location and element of each atom in the structure. For bulk structures, we typically get the initial structure file from <a href="https://materialsproject.org/" target="_blank">The Materials Project</a>. My preferred method is to just Google for my desired material (e.g. “InAs Materials Project”). Once you have the bulk structure, slabs or interfaces can be generated using either VaspVis or Ogre (our groups Python packages). More details about this with code examples will come later. The only calculations that require additional alterations to the POSCAR file are for an unfolded Band calculation and an OPT calculation

### POTCAR
The POTCAR is dependent on the elements used in the calculation and the order of the elements in the POSCAR. It is crucial that the order of the elements in the POTCAR match the order of the elements in the sixth line of the POSCAR. A POTCAR can be easily generated using the potcar.sh file. For example a common way to generate the POTCAR file is to use the head command to view the top lines of the POSCAR, and then use the potcar.sh file to generate the POTCAR file using the same elements shown in the 6th line of the POSCAR. To check the POTCAR you can use the grep command.

Given a POSCAR you can create the correct POTCAR using the following commands:
```
(base) ddardzin(cori)$ head POSCAR
In1 As1
1.0
0.000000 3.029200 3.029200
3.029200 0.000000 3.029200
3.029200 3.029200 0.000000
In As
1 1
direct
0.000000 0.000000 0.000000 In
0.750000 0.750000 0.750000 As
(base) ddardzin(cori)$ potcar.sh In As
(base) ddardzin(cori)$ grep 'TITEL' POTCAR
	TITEL = PAW_PBE In 08Apr2002  
	TITEL = PAW_PBE As 22Sep2009
```

## Self-Consistent Field (SCF) Calculation
The SCF calculation will be used to generate converged CHG and CHGCAR files so they can be used to calculate the eigenvalues of the system.

### INCAR
In the INCAR there are four important parameters beyond the general parameters. The one parameter that never changes is ICHARG=2 because the SCF calculation must calculate the charge density of the system. The other three parameters can change based on your needs. Sometimes we don’t need to write the CHG* files since we only need to do the SCF calculation, or sometimes we need the WAVECAR for STM simulations or wave function visualizations. Additionally, ISMEAR=0 is the standard parameter we used, but for certain materials a different smearing method might be necessary for convergence.

```txt
ICHARG = 2 # Generate CHG* from a superposition of atomic charge densities
ISMEAR = 0 # Fermi smearing
LCHARG = True # Write the CHG* files
LWAVE = False # Does not write the WAVECAR
```

The following code can be used to generate the INCAR
```bash
incar.py --scf
```
or
```bash
incar.py -s
```

### KPOINTS
The KPOINTS file for an SCF calculation is one of the most basic and easy to generate  files. For the SCF calculation we use a Γ-centered Monkhorst-Pack grid of a specified density. Depending on the material, the size of the grid might change. For example, if you are working with a supercell, the k-point grid does not have to be as dense since the reciprocal lattice size is inversely proportional to the real space lattice size (i.e. large real space lattice $\rightarrow$ small reciprocal lattice $\rightarrow$ less dense grid required to fill the space)

```txt
Automatic kpoint scheme  
0  
Gamma  
7 7 7
```

To generate the KPOINTS file for a SCF calculation the kpoints.py file can be used.
```bash
kpoints.py --grid --density 7 7 7
```
or
```bash
kpoints.py -g -d 7 7 7
```

## Density of States (DOS) Calculation
### INCAR
In the INCAR there are eight important parameters beyond the general parameters. The only parameters that can be changed are NEDOS, EMIN, and EMAX. These values are used to determine the energy range to be calculated and how many points are included within the given energy range. By default, only 301 points are sampled for the density of states, and the energy range is unconstrained, so it includes the states with the lowest and highest energies. For the most part, we only care about a small range of energies around the Fermi energy, so we like to constrain the values between emin and emax.