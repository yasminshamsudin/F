# 2nd GEN electric fields calculation scripts

The Second Generation scripts use a pdb file as starting structure and generate AM1-BCC charges from antechamber and creates gromacs-compatible parameter and structure files.

In 1.pdb2amber2gmx.sh a pdb file gets GAFF vdw parameters and electrostatic parameters from AM1-BCC. 
It is then converted to gromacs .top and .gro files using the acpype.py script.

In 2a.solvateGAFF.sh the user specifies the solvent box. It's possible to specify a range of sizes for testing.
The script generates the box and fills it with solvents specified by the user based on the densities provided in the densities.F file.
The 2b.solvateOPLS.sh is similar to the GAFF counterpart, but adds one more line in the header to include the opls forcefield.itp file defining the OPLS forcefield.

In 3.createScripts.sh all the input files for simulations in GROMACS and the analysis script are generated.
The 3b.createComplexScripts.sh is used for systems with more than one ligand (like a ligand and a counterion).

The 4.runGMXsherlock.sh should be placed on the Sherlock cluster (or modified to suit whichever cluster the user is familiar with).
It submits the simulations to Sherlock using sbatch.

The 5.get.Fields.sh summarizes the data from all simulations into one log-file for easy analysis of electric fields, standard deviations, standard error of the means, etc.
