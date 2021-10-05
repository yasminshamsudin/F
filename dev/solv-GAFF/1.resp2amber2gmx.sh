#!/bin/bash

# By Yasmin 2021-06-24 
# 2021-07-07 Tested locally 
# 2021-07-20 Integrated into F
# 2021-07-20 Added functionality to check for charged ligands and add counterions if charged
# 2021-07-20 Tested locally
# 2021-08-11 Formatting updates
# 2021-08-18 Formatting updates

# This script creates ligand AMBER format toppology and lib files. 
# AMBER top and lib files are then converted to GROMACS gro and top files.

# Requires an AMBER antechamber installation to get AMBER files.
# Requires acpype.py for conversion to GROMACS formatted files.

####################### EDIT ONLY THESE PARAMETERS #########################################

molType=LIG                             # Choose SOL for solvents, LIG for ligands

runDirectory=RESP                            # Path to the run directory
pythonVersion=python3.8                       # Installed python version

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

AMBERPATH=/home/yasmin/amber20/amber.sh     # Path for the Amber installation
FFPATH=/home/yasmin/amber20_src/dat/leap/parm/gaff2.dat   # default: $AMBERHOME/dat/leap/parm/gaff2.dat
ligFF=leaprc.gaff2

######################## NO MORE EDITING PAST THIS LINE! ###################################


# Set the AMBER path to find all modules

source $AMBERPATH

# Convert ligand pdb to mol2 using antechamber

cd $runDirectory

for i in *.resp
 do
    ligname=$( echo "$i" | sed -e 's/_opt\.resp//g')

    # Run antechamber to generate .in-file
        antechamber -i $i -fi gesp -o $ligname.in -fo prepi -c resp -rn $molType -pf y 
        wait
done

# Check the parameters using parmchk2

for i in *.in
 do
    ligname=$( echo "$i" | sed -e 's/\.in//g')
    parmchk2 -i $ligname.in -f prepi -o $ligname.frcmod -a Y -p $FFPATH
    wait
done

# Create the .lib and .prmtop files using tleap

for i in *.frcmod
 do
    ligname=$( echo "$i" | sed -e 's/\.frcmod//g')
     echo -e "source $ligFF" > tleap.$ligname.in
     echo -e "loadAmberprep $ligname.in" >> tleap.$ligname.in
     echo -e "loadamberparams $ligname.frcmod" >> tleap.$ligname.in
     echo -e "check $molType" >> tleap.$ligname.in
     echo -e "saveoff $molType $ligname.lib" >> tleap.$ligname.in
     echo -e "saveamberparm $molType $ligname.prmtop $ligname.inpcrd" >> tleap.$ligname.in
     echo -e "Quit" >> tleap.$ligname.in

    tleap -f tleap.$ligname.in
    wait
done

# Convert to GROMACS files using acpype.py (uses GAFF parameters and ff99SB forcefield)

for i in *.prmtop
 do
   ligname=$( echo "$i" | sed -e 's/\.prmtop//g')
   $pythonVersion $workingDirectory/acpype.py -p $workingDirectory/$runDirectory/$ligname.prmtop -x $workingDirectory/$runDirectory/$ligname.inpcrd
   wait

   # Clean up folder & filenames and ligand names in files
   mkdir $ligname
   mv $molType"_GMX.top" $ligname/$ligname"_GMX.top"
   mv $molType"_GMX.gro" $ligname/$ligname"_GMX.gro"
   
   mkdir $ligname/AMBER
   mv *$ligname.* $ligname/AMBER
done

rm *.*
 
###### DONE! ##########
