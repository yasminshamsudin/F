#!/bin/bash

# By Yasmin 2021-07-02 
# 2021-07-12 Added replicates functionality 
# 2021-08-11 Integrated into F
# 2021-08-11 Added functionality for adaptable computational time based on box size
# 2021-08-11 Tested locally
# 2021-08-12 Tested on Sherlock
# 2021-08-18 Formatting changes

# This script prepares job files for GROMACS simulations. 
# This script submits files on the Stanford Sherlock cluster.
# The submit file also runs a script for calculating electric fields. 

# Requires access on Stanford's Sherlock cluster. Can be modified for other clusters.

############################# EDIT THESE PARAMETERS #########################################

runDirectory=211005-SAM             # Path to the ligand directory
pythonVersion=python                        # Installed python version

partitions=hns,iric,owners                  # Default for Boxers: hns

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

nodes=1                                     # Default: 1
cores=8                                     # Default: 8 on Sherlock
replicates=5				    # Default:1

######################## NO MORE EDITING PAST THIS LINE! ###################################

cd $runDirectory

# Repeat for all ligands
for ligname in *
 do
    cd $ligname

    # Create generic job-files

        # Assign Sherlock flags for running
        echo -e "#!/bin/bash" > COMPLEX.job
        echo -e "#SBATCH --mail-type=END" >> COMPLEX.job
        echo -e "#SBATCH --job-name=COMPLEX" >> COMPLEX.job
        echo -e "#SBATCH -p $partitions" >> COMPLEX.job
        echo -e "#SBATCH -N $nodes" >> COMPLEX.job
        echo -e "#SBATCH --ntasks-per-node=1" >> COMPLEX.job
        echo -e "#SBATCH --cpus-per-task=$cores" >> COMPLEX.job
        echo -e "#SBATCH --time=WALLTIME"'\n' >> COMPLEX.job

        # Load latest version of gromacs and python 2.7 modules
        echo -e "module reset" >> COMPLEX.job
	    echo -e "module load chemistry gromacs" >> COMPLEX.job
        echo -e "module load py-numpy/1.14.3_py27"'\n' >> COMPLEX.job

        # Run minimization
        echo -e "gmx grompp -f ../min.mdp -c COMPLEX.gro -p COMPLEX.top -o COMPLEX_min.tpr" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_min"'\n' >> COMPLEX.job

        # Run heating
        echo -e "gmx grompp -f ../nvt.mdp -c COMPLEX_min.gro -p COMPLEX.top -o COMPLEX_nvt.tpr -maxwarn 2" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_nvt"'\n' >> COMPLEX.job

        # Run equilibration
        echo -e "gmx grompp -f ../npt.mdp -c COMPLEX_nvt.gro -p COMPLEX.top -o COMPLEX_npt.tpr -maxwarn 2" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_npt"'\n' >> COMPLEX.job

        # Run data collection
        echo -e "gmx grompp -f ../md.mdp -c COMPLEX_npt.gro -p COMPLEX.top -o COMPLEX_md.tpr -maxwarn 2" >> COMPLEX.job
        echo -e "gmx mdrun -v -deffnm COMPLEX_md"'\n' >> COMPLEX.job

        # Extract coordinates and forces from simulation
        echo -e "gmx traj -s COMPLEX_md.tpr -f COMPLEX_md.trr -n ../probe.ndx -ox co_coords.xvg" >> COMPLEX.job
        echo -e "gmx traj -s COMPLEX_md.tpr -f COMPLEX_md.trr -n ../probe.ndx -of co_forces.xvg" >> COMPLEX.job

        # Re-run data collection with solvent charges neutralized
        echo -e "gmx grompp -f ../md.mdp -c COMPLEX_npt.gro -p COMPLEX_0q.top -o COMPLEX_md_0q.tpr -maxwarn 2" >> COMPLEX.job
        echo -e "gmx mdrun -rerun COMPLEX_md.trr -v -deffnm COMPLEX_md_0q"'\n' >> COMPLEX.job

        # Extract coordinates and forces from simulation of uncharged environment
        echo -e "gmx traj -s COMPLEX_md_0q.tpr -f COMPLEX_md_0q.trr -n ../probe.ndx -ox co_coords_0q.xvg" >> COMPLEX.job
        echo -e "gmx traj -s COMPLEX_md_0q.tpr -f COMPLEX_md_0q.trr -n ../probe.ndx -of co_forces_0q.xvg"'\n' >> COMPLEX.job
        wait

        # Run analysis script to extract electric fields
        echo -e "cp ../BLCO.py ." >> COMPLEX.job
        echo -e "python BLCO.py" >> COMPLEX.job
        echo -e "wait"'\n' >> COMPLEX.job
        echo -e "rm BLCO.py" >> COMPLEX.job

    # Customize job-files for all solvents

    for i in *0q.top   # Get the name of each solvent
        do

	   # Create replicates
	   for (( n = 1; n <= $replicates; n++ ))
	     do
                complex=$( echo "$i" | sed -e 's/\_0q.top//g')
                mkdir $complex"_s"$n
                sed "s/COMPLEX/$complex/g" COMPLEX.job > $complex"_s"$n/run.job
                cp $complex.* $complex"_s"$n/	
                cp $complex*.top $complex"_s"$n/
                cd $complex"_s"$n

                # Check the box size
                getboxsize=$(echo "$complex" | tr -dc '0-9')
                boxsize=$(echo "$getboxsize" | bc)

               # If box size is larger than 20:
               if [ $boxsize -gt 20 ] && [ $boxsize -le 30 ]
                       then
                       walltime=24:00:00 # most simulations finish in 3 h
               elif [ $boxsize -gt 30 ] && [ $boxsize -le 40 ]
                       then
                       walltime=48:00:00 # most simulations finish in 5 h
               elif [ $boxsize -gt 40 ]
                       then
                       walltime=48:00:00
               # Otherwise for boxsizes < 20:
               else
                       walltime=05:00:00 # most simulations finish in 1 h
               fi

                sed -i "s/WALLTIME/$walltime/g" run.job

                # Submit the run on Sherlock
	        sbatch run.job
                cd ..

	   done # End for all replicates

           # Housekeeping
           rm $complex.*
	   rm $complex*.top

        done # End for all solvents
  
        cd ..

done # End repeat for all ligands
