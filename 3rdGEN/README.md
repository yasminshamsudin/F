# 3rd GEN scripts for electric field calculations

3rd GEN scripts for preparing and running electric field calculations.
These scripts prepare ligand files starting with RESP files from GAUSSIAN opt/freq + resp simulations.
1. RESP files are converted to AMBER format through antechamber and then to GROMACS format through acpype.py
2. Ligands are solvated in various solvents and in user-defined box sizes.
3. Input files for simulations and analyses are created.
4. Files are submitted for simulations on clusters.
