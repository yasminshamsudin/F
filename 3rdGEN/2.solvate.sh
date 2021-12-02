#!/bin/bash

# By Yasmin 2021-09-30 
# 2021-09-30 Adapted from script in Github repository (2nd Gen)
# 2021-12-02 Bug fixes

# This script solvates ligands in a chosen solvent for GROMACS simulations.

# Requires GROMACS installation.

############################# EDIT THESE PARAMETERS #########################################

FF=o5r                     # Forcefield options: g (GAFF, default), o (OPLSAA), a (AMBER 14SB)
runDirectory=PREP          # Path to the ligand directory

############ Advanced modelling. Only edit if you know what you are doing! #################

workingDirectory=$(pwd)         # Path to where this script is located
solventDirectory=SOLVENTS       # Path to directory of solvents

boxtype=cubic                   # default solvent box: cubic
boxradiusSmallest=20            # default radius (in Å): 20
boxradiusLargest=40
boxsizeSteps=10

water=tip3p				        # default water model: tip3p

######################## NO MORE EDITING PAST THIS LINE! ###################################

# Go to ligand directory
cd $runDirectory

# Repeat for all ligands
for ligname in *
 do
    # Create log file to keep track of solvation
    echo -e "Solvent solvated expected percentage" > ../$ligname.solvate.log
    
	cd $ligname

    # Repeat for all box sizes
    for ((boxsize = $boxradiusSmallest; boxsize <= $boxradiusLargest; boxsize=$(( $boxsize + $boxsizeSteps ))))
     do

        # Convert the radius to nm for gmx:
        gmxsize=$(( $boxsize / 10 ))

        # Center the solute in the box
        gmx editconf -f $ligname"_GMX.gro" -o $ligname"_"$boxsize.gro -bt $boxtype -d $gmxsize
        wait

        # Check the number of molecules in the file
        awk '{print $1}' $ligname"_GMX.gro" > noligands.tmp
        tail -2 noligands.tmp > twonoligands.tmp
        noligands=$( head -1 twonoligands.tmp) 
        
        # Create ligand in other solvents (non-water) boxes
        
        cd $workingDirectory/$solventDirectory/

        for i in *.gro   # Repeat for all solvents in the SOLVENTS folder
          do
            solvname=$( echo "$i" | sed -e 's/_GMX\.gro//g')
            cd $workingDirectory/$runDirectory/$ligname/

            # Calculate number of molecules needed in the box to get correct density

                # Get the parameters for the solvent
                grep "$solvname" $workingDirectory/densities.F > $solvname"_density.tmp"

                # Get MW from the 2nd column and scale to avoid decimals
                awk '{print $2 * 1000}' $solvname"_density.tmp" > MW.tmp
                scaledMW=$( tail MW.tmp)

                # Get the density from the 3rd column and scale to avoid decimals
                awk '{print $3 * 1000}' $solvname"_density.tmp" > density.tmp
                scaledDensity=$( tail density.tmp)

                # Use avogadro scaled to avoid decimals
                scaledAvogadro=6022

                # Clean up the folder
                rm *.tmp

                # Calculate the number of solvent molecules needed to fill the box
                solventMolecules=$(($scaledDensity * $boxsize * 2 * $boxsize * 2 * $boxsize * 2 * $scaledAvogadro / $scaledMW))
                scaledSolvents=$(($solventMolecules / 10000))
                nsolvmol=$(echo ${scaledSolvents%.*})
                echo "number of solvent molecules: $nsolvmol"

	        # Solvate the box with solvent
            gmx insert-molecules -f $ligname"_"$boxsize.gro -ci $workingDirectory/$solventDirectory/$i -nmol $nsolvmol -o ${ligname}_${solvname}${FF}_${boxsize}.gro
            wait
                
                # Check that the box has the correct number of solvents
                awk '{print $1}' ${ligname}_${solvname}${FF}_${boxsize}.gro > noligands.tmp
                tail -2 noligands.tmp > twonoligands.tmp
                getnoligands=$(head -1 twonoligands.tmp )
                nosystem=$(echo "$getnoligands" | tr -dc '0-9')
                addedsolvents=$(( $nosystem - $noligands ))

                # If the box is not properly solvated, try again
                    if [ $addedsolvents -ne $nsolvmol ]
                            then
                            gmx insert-molecules -f $ligname"_"$boxsize.gro -ci $workingDirectory/$solventDirectory/$i -nmol $nsolvmol -rot none -o ${ligname}_${solvname}${FF}_${boxsize}.gro
                            rm *.gro.*

                            awk '{print $1}' ${ligname}_${solvname}${FF}_${boxsize}.gro > noligands.tmp
                            tail -2 noligands.tmp > twonoligands.tmp
                            getnoligands=$(head -1 twonoligands.tmp )
                            nosystem=$(echo "$getnoligands" | tr -dc '0-9')
                            addedsolvents=$(( $nosystem - $noligands))

                    fi

                # Print the number of solvents and the percentage of the desired number of solvents in the box
                percentage=$(( $addedsolvents * 100 / $nsolvmol))
                echo -e "$ligname.$solvname $addedsolvents $nsolvmol $percentage " >> $workingDirectory/$ligname.solvate.log

	        # Create master topology (.top) files 

                # Get the header until atom types from ligand file
                sed -n -r '1,/atomtypes/p' $ligname"_GMX.top" >> header.tmp    # Copy all lines until atomtypes

                # Get the atomtypes of ligand and solvent and remove duplicates
                cp $workingDirectory/$solventDirectory/$solvname$FF.itp .    # Copy solvent file to ligand directory
                sed -n -r '/atomtypes/,/^\s*$/p' $solvname$FF.itp > solvatomtypes.tmp    # Get atomtypes of solvent
                sed -n -r '/atomtypes/,/^\s*$/p' $ligname"_GMX.top" > ligatomtypes.tmp    # Get atomtypes of ligand
                sed -i '1,2d;$d' *atomtypes.tmp                                     # Remove first and last (empty) lines
                cat solvatomtypes.tmp ligatomtypes.tmp > complexatomtypes.tmp    # Put it together
                sort complexatomtypes.tmp | uniq > complex_sorted.tmp                   # Remove duplicate lines
                echo -e "" >> complex_sorted.tmp                                   # Adds blank line at the end

                # Get ligand parameters and include solvent parameters
                sed -n -r '/moleculetype/,/system/p' $ligname"_GMX.top" > ligparams.tmp # Get the ligand parameters
                sed -i '$ d' ligparams.tmp                                     # Remove first and last (empty) lines
#                echo -e "#include \"../ions.itp\"\n" >> ligparams.tmp # Add ion parameters
                echo -e "#include \"../$solvname$FF.itp\"\n" >> ligparams.tmp # Add solvent parameters

                sed -n -r '/system/,//p' $ligname"_GMX.top" > system.tmp # Get the system name parameters
                echo -e " SOL              $addedsolvents" >> system.tmp # Add the number of solvent molecules

                # Put everything together into top files (charged and no-charge)
                cat header.tmp complex_sorted.tmp ligparams.tmp system.tmp > ${ligname}_${solvname}${FF}_${boxsize}.top # Create charged top-file
#                sed -i 's/MOL/LIG/g' ${ligname}_${solvents}_${boxsize}.top     # Change all instances of MOL to LIG

                cp ${ligname}_${solvname}${FF}_${boxsize}.top ${ligname}_${solvname}${FF}_${boxsize}"_0q.top" # Make no-charge version
                sed -i "s/$solvname$FF/$solvname$FF"_0q"/g" ${ligname}_${solvname}${FF}_${boxsize}"_0q.top"   # Change path to no-charge solvent itp  
            
	        # Create parameter files for charged and uncharged solvents (.itp)
            
                # Create a solvent itp file from file copied from source for importing
                sed -i '/moleculetype/,$!d' $solvname$FF.itp                   # Make the charged solvent itp file

                # Create a no-charge version for the atom types
                sed -n -r '/moleculetype/,/mass/p' $solvname$FF.itp > header_0q.tmp
                sed -n -r '/mass/,/^\s*$/p' $solvname$FF.itp > solvatoms.tmp    # Get atoms section of solvent
                sed -i '1d;$d' solvatoms.tmp                                     # Remove first line and blank last line
                awk '{print "     "$1"         "$2"      "$3"    "$4"     "$5"      "$6"    0.00000  "$8}' solvatoms.tmp  > solvatoms_0q.tmp
                echo -e "" >> solvatoms_0q.tmp                                   # Adds blank line at the end

                sed -n -r '/bonds/,//p' $solvname$FF.itp > end.tmp         # Get the system name parameters                      
                cat header_0q.tmp solvatoms_0q.tmp end.tmp > $solvname$FF"_0q.itp"    # Create no-charge solvent itp file

	        # Clean up the folder

            rm *.tmp                                            
            wait
            cd $workingDirectory/$solventDirectory/

        done

        # Set up water calculations for each box size

            cd $workingDirectory/$runDirectory/$ligname/

            # Solvate in water	
	        cp $ligname"_GMX.top" ${ligname}_${water}.tmp.top
	        gmx solvate -cp $ligname"_"$boxsize.gro -cs spc216.gro -o ${ligname}_${water}_${boxsize}.gro -p ${ligname}_${water}.tmp.top

	        # Add the forcefield and water model itp file to topology (.top) file.	
            echo -e "#include \"../forcefield.itp\"\n" > ${ligname}_${water}_${boxsize}.top # Add forcefield parameters	
            sed -n -r '/atomtypes/,/system/p' ${ligname}_${water}.tmp.top >> ${ligname}_${water}_${boxsize}.top # Add the main body
            sed -i '$d' ${ligname}_${water}_${boxsize}.top # Removes the last line
#            echo -e "#include \"../ions.itp\"\n" >> ${ligname}_${water}_${boxsize}.top # Add ion parameters
            echo -e "#include \"../$water.itp\"\n" >> ${ligname}_${water}_${boxsize}.top # Add solvent parameters
            sed -n -r '/system/,//p' ${ligname}_${water}.tmp.top >> ${ligname}_${water}_${boxsize}.top # Get the system name parameters until end of file
            #sed -i 's/MOL/LIG/g' ${ligname}_${water}_${boxsize}.top     # Change all instances of MOL to LIG
            rm *.tmp.*
    
            # Copy the right solvent files to the directory
            cp $workingDirectory/oplsaa.ff/$water.itp .
            cp $workingDirectory/oplsaa.ff/forcefield.itp .
            cp $workingDirectory/oplsaa.ff/ff*.itp .
            cp $workingDirectory/oplsaa.ff/ions.itp .

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

    done # End repeat for all box sizes
    
    # Clean up
    mkdir GMXPREP
    mv $ligname"_GMX".* GMXPREP/
    rm *.tmp
    cd $workingDirectory/$runDirectory

done # End repeat for all ligands
