# Calculation Descriptions
In this section we will discus the inputs you need in order to run all the different types of calculations that we run in the Marom group.

## General Inputs
### INCAR
As stated in [[VASP Basics]], the INCAR contains the set of instructions that tells VASP what type of  calculation to perform. Although we perform many different types of calculations, there are some standard parameters that we always use. Sometimes, we may choose to alter   the following parameters slightly, but for the vast majority of our calculations, these will remain constant.

> ```bash
	ALGO = Fast # Mixture of Davidson and RMM-DIIS algos  
	PREC = N # Normal precision  
	EDIFF = 1e-5 # Convergence criteria for electronic converge  
	NELM = 500 # Max number of electronic steps  
	ENCUT = 350 # Cut off energy  
	LASPH = True # Include non-spherical contributions from gradient corrections  
	NBANDS = 100 # Number of bands to include in the calculation  
	BMIX = 3 # Mixing parameter for convergence  
	AMIN = 0.01 # Mixing parameter for convergence  
	SIGMA = 0.05 # Width of smearing in eV
```


