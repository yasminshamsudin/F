#!/bin/bash

# By Yasmin 2021-06-24 
# 2021-07-07 Tested locally
# 2021-07-09 Tested on Sherlock
# 2021-07-13 Added water setup
# 2021-07-24 Integrated into F
# 2021-08-05 Changed from gmx solvate to insert-molecules to add boxsize flexibility
# 2021-08-06 Tested locally & verified correct densities

# This script solvates ligands in a chosen solvent for GROMACS simulations.

# Requires GROMACS installation.

############################# EDIT THESE PARAMETERS #########################################

FF=GAFF-Singles                                # default force field: GAFF-Yasmin

ligandDirectory=RUN                     # Path to the ligand directory
solventDirectory=SOLVENTS                   # Path to directory of solvents

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)                     # Path to where this script is located

boxtype=cubic                               # default solvent box: cubic
boxsizeSmallest=30                         # default radius (in Å): 20
boxsizeLargest=40
boxsizeSteps=10
noligands=1                                 # default number of ligands: 1

water=tip3p				                    # default water model: tip3p

######################## NO MORE EDITING PAST THIS LINE! ###################################

# Activate GROMACS libraries
#ml reset
#ml load chemistry gromacs

# Go to ligand directory
cd $ligandDirectory

# Repeat for all ligands
for ligname in *
 do
	cd $ligname

    # Repeat for all box sizes
    for ((boxsize = $boxsizeSmallest; boxsize <= $boxsizeLargest; boxsize=$(( $boxsize + $boxsizeSteps ))))
     do

        # Convert the radius to nm for gmx:
        gmxsize=$(( $boxsize / 10 ))

        # Center the solute in the box
        gmx editconf -f $ligname"_GMX.gro" -o $ligname"_"$boxsize.gro -bt $boxtype -d $gmxsize
        wait

        # Create ligand in other solvents (non-water) boxes
        # Repeat for all solvents in the SOLVENTS folder

        cd $workingDirectory/$solventDirectory/$FF/

        for solvents in *   # Get the name of each solvent
         do
            cd $solvents
                for i in *.pdb
                do
                solvname=$( echo "$i" | sed -e 's/\.pdb//g')
                cd $workingDirectory/$ligandDirectory/$ligname/

                # Calculate number of molecules needed in the box to get correct density

                    # Get the parameters for the solvent
                    grep "$solvents" $workingDirectory/densities.txt > $solvents"_density.tmp"

                        # Get MW from the 2nd column and scale to avoid decimals
                        awk '{print $2 * 1000}' $solvents"_density.tmp" > MW.tmp
                        scaledMW=$( tail MW.tmp)

                        # Get the density from the 3rd column and scale to avoid decimals
                        awk '{print $3 * 1000}' $solvents"_density.tmp" > density.tmp
                        scaledDensity=$( tail density.tmp)

                        # Use avogadro scaled to avoid decimals
                        scaledAvogadro=6022

                        # Clean up the folder
                        rm *.tmp

                    # Calculate the number of solvent molecules needed to fill the box
                    solventMolecules=$(($scaledDensity * $boxsize * $boxsize * $boxsize * $scaledAvogadro /$scaledMW))
                    scaledSolvents=$(($solventMolecules / 10000))
                    nsolvmol=$(echo ${scaledSolvents%.*})
                    echo "number of solvent molecules: $nsolvmol"

	            # Solvate the box with solvent

                    gmx insert-molecules -f $ligname"_"$boxsize.gro -ci $workingDirectory/$solventDirectory/$FF/$solvents/$solvname.pdb -nmol $nsolvmol -o ${ligname}_${solvents}_${boxsize}.gro
                    wait
                

	            # Create master topology (.top) files 

                # Get the header until atom types
                sed -n -r '/MOL_GMX/,/Amb/p' $ligname"_GMX.top" > header.tmp    

                # Get the atomtypes of ligand and solvent and remove duplicates
                cp $workingDirectory/$solventDirectory/$FF/$solvents/$solvname.itp .    # Copy solvent file to ligand directory
                sed -n -r '/atomtypes/,/^\s*$/p' $solvname.itp > solvatomtypes.tmp    # Get atomtypes of solvent
                sed -n -r '/atomtypes/,/^\s*$/p' $ligname"_GMX.top" > ligatomtypes.tmp    # Get atomtypes of ligand
                sed -i '1,2d;$d' *atomtypes.tmp                                     # Remove first and last (empty) lines
                cat solvatomtypes.tmp ligatomtypes.tmp > complexatomtypes.tmp    # Put it together
                sort complexatomtypes.tmp | uniq > complex_sorted.tmp                   # Remove duplicate lines
                echo -e "" >> complex_sorted.tmp                                   # Adds blank line at the end

                # Get ligand parameters and include solvent parameters
                sed -n -r '/moleculetype/,/system/p' $ligname"_GMX.top" > ligparams.tmp # Get the ligand parameters
                sed -i '$ d' ligparams.tmp                                     # Remove first and last (empty) lines
                echo -e "#include \"../$solvname.itp\"\n" >> ligparams.tmp # Add solvent parameters

                sed -n -r '/system/,//p' $ligname"_GMX.top" > system.tmp # Get the system name parameters
                echo -e " SOL              $nsolvmol" >> system.tmp # Add the number of solvent molecules

                # Put everything together into top files (charged and no-charge)
                cat header.tmp complex_sorted.tmp ligparams.tmp system.tmp > ${ligname}_${solvents}_${boxsize}.top # Create charged top-file
                sed -i 's/MOL/LIG/g' ${ligname}_${solvents}_${boxsize}.top     # Change all instances of MOL to LIG

                cp ${ligname}_${solvents}_${boxsize}.top ${ligname}_${solvents}_${boxsize}"_0q.top" # Make no-charge version
                sed -i "s/$solvname/$solvname"_0q"/g" ${ligname}_${solvents}_${boxsize}"_0q.top"   # Change path to no-charge solvent itp  
            
	            # Create parameter files for charged and uncharged solvents (.itp)
            
                # Create a solvent itp file from file copied from source for importing
                sed -i '/moleculetype/,$!d' $solvname.itp                   # Make the charged solvent itp file

                # Create a no-charge version for the atom types
                sed -n -r '/moleculetype/,/mass/p' $solvname.itp > header_0q.tmp
                sed -n -r '/mass/,/^\s*$/p' $solvname.itp > solvatoms.tmp    # Get atoms section of solvent
                sed -i '1d;$d' solvatoms.tmp                                     # Remove first line and blank last line
                awk '{print "     "$1"         "$2"      "$3"    "$4"     "$5"      "$6"    0.00000  "$8}' solvatoms.tmp  > solvatoms_0q.tmp
                echo -e "" >> solvatoms_0q.tmp                                   # Adds blank line at the end

                sed -n -r '/bonds/,//p' $solvname.itp > end.tmp         # Get the system name parameters                      
                cat header_0q.tmp solvatoms_0q.tmp end.tmp > $solvname"_0q.itp"    # Create no-charge solvent itp file

	            # Clean up the folder

                rm *.tmp                                            
                wait
                cd $workingDirectory/$solventDirectory/$FF/
                done
        done

        # Set up water calculations for each box size

        cd $workingDirectory/$ligandDirectory/$ligname/

        # Solvate in water
	
	    cp $ligname"_GMX.top" ${ligname}_${water}.tmp.top
	    gmx solvate -cp $ligname"_"$boxsize.gro -cs spc216.gro -o ${ligname}_${water}_${boxsize}.gro -p ${ligname}_${water}.tmp.top

	    # Add the forcefield and water model itp file to topology (.top) file.
 	
       echo -e "#include \"../forcefield.itp\"\n" > ${ligname}_${water}_${boxsize}.top # Add forcefield parameters	
       sed -n -r '/atomtypes/,/system/p' ${ligname}_${water}.tmp.top >> ${ligname}_${water}_${boxsize}.top # Add the main body
       sed -i '$d' ${ligname}_${water}_${boxsize}.top # Removes the last line
       echo -e "#include \"../$water.itp\"\n" >> ${ligname}_${water}_${boxsize}.top # Add solvent parameters
        sed -n -r '/system/,//p' ${ligname}_${water}.tmp.top >> ${ligname}_${water}_${boxsize}.top # Get the system name parameters until end of file
        sed -i 's/MOL/LIG/g' ${ligname}_${water}_${boxsize}.top     # Change all instances of MOL to LIG
        rm *.tmp.*
    
        # Copy the right solvent files to the directory

        cp $workingDirectory/oplsaa.ff/$water.itp .
        cp $workingDirectory/oplsaa.ff/forcefield.itp .
        cp $workingDirectory/oplsaa.ff/ff*.itp .

        # Create a no-charge versions for the water runs
            cp ${ligname}_${water}_${boxsize}.top ${ligname}_${water}_${boxsize}"_0q.top" # Make no-charge version
            sed -i "s/$water/$water"_0q"/g" ${ligname}_${water}_${boxsize}"_0q.top"   # Change path to no-charge solvent itp

            sed -n -r '/moleculetype/,/charge/p' $water.itp > header_0q.tmp
            sed -n -r '/charge/,/^\s*$/p' $water.itp > solvatoms.tmp    # Get atoms section of solvent
            sed -i '1d;$d' solvatoms.tmp                                     # Remove first line and blank last line
            awk '{print $1"  "$2"      "$3"       "$4"             "$5"             "$6"        0.00000  "}' solvatoms.tmp  > solvatoms_0q.tmp
            echo -e "" >> solvatoms_0q.tmp                                   # Adds blank line at the end

            sed -n -r '/FLEXIBLE/,//p' $water.itp > end.tmp         # Get the system name parameters                      
            cat header_0q.tmp solvatoms_0q.tmp end.tmp > $water"_0q.itp"    # Create no-charge solvent itp file

    done
    
        # Clean up

         mkdir GMXPREP
         mv $ligname"_GMX".* GMXPREP/
         mv $ligname.gro GMXPREP/
         cd $workingDirectory/$ligandDirectory
done