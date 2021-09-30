#!/bin/bash

# By Yasmin 2021-07-02 
# 2021-07-02 Tested locally 
# 2021-07-05 Tested on Sherlock 
# 2021-07-14 Added capability to process several replicates
# 2021-07-14 Added capability to merge several ligand runs into one file
# 2021-07-14 Tested on Sherlock 
# 2021-08-16 Added loop for iterating over box sizes
# 2021-08-19 Tested locally
# 2021-08-26 Tested on Sherlock
# 2021-09-30 Stopped development of script and deposited to GitHub.

# This script analyzes electric field calculations from GROMACS simulations. 

############################# EDIT THESE PARAMETERS #########################################

runNAME=avboh15                  # Your identifier for this simulation
runDirectory=RUN15                     # Path to the simulations directory

############ Advanced modelling. Only edit if you know what you are doing! #################

date=$( date +%y%m%d)                       # Today's date
identifier=$date-$runNAME                 # Unique identifiers for the group of simulations
workingDirectory=$(pwd)                     # Path to where this script is located

######################## NO MORE EDITING PAST THIS LINE! ###################################

echo -e "Results of $identifier" > $identifier.header.tmp
cp $identifier.header.tmp summary.tmp

cd $runDirectory

# Repeat for all ligands
for ligname in *
 do
    cd $ligname

    # Repeat for all box sizes
    for getboxsize in ${ligname}_*.gro
        do
        boxsizename=$(echo "$getboxsize" | sed -e "s/${ligname}_//g")
        boxsize=$(echo "$boxsizename" | tr -dc '0-9')
        echo -e "Boxradius of $boxsize" >> $workingDirectory/summary.tmp

        # Create titles for each column in file for each replicate
        echo -e "Solvent <Field> Standard_deviation Removed_lines" > $ligname.header.tmp
        cp $ligname.header.tmp combinedfields.tmp

        # Create titles for each column in file combining all replicates
        echo -e "${ligname}_Solvent <Field> Standard_error <StdDev> Data_points" > $ligname.header.tmp
        cp $ligname.header.tmp final.tmp

        # Repeat for all solvents
        for solvents in *_0q.itp
	      do
           solvent=$( echo "$solvents" | sed -e 's/_0q\.itp//g')
           
           # Initiate for calculating averages and SEM over replicates
            echo -e "" > $solvent.replicate.tmp
            cp $solvent.replicate.tmp $solvent.stdDev.replicate.tmp
            cp $solvent.replicate.tmp $solvent.counted.tmp
            
           # Initiate by making the first column with the solvent name
           echo -e "$solvent" > $ligname.$solvent.name.tmp
           cp $ligname.$solvent.name.tmp combinedsolvents.tmp

            # Repeat for all replicates of the solvent
            for solvreplica in ${ligname}_${solvent}_${boxsize}_s*
             do
                cd $solvreplica

                # Remove lines where probe is at the edge of the box and count the number of removed lines
                awk '$6 <= 1' BLCO_md.txt > counted_lines.txt
                awk '$6 > 1' BLCO_md.txt > removed_lines.txt
                awk 'END{print NR}' removed_lines.txt > ../$ligname.$solvent.noremovedlines.tmp
                awk 'END{print NR}' counted_lines.txt > ../$ligname.$solvent.nocountedlines.tmp
                head -1 ../$ligname.$solvent.nocountedlines.tmp >> ../$solvent.counted.tmp # Output the value to calculate average of averages

                # Calculate the average fields
                awk '{ total += $4 } END { print total/NR }' counted_lines.txt > ../$ligname.$solvent.aveField.tmp
                head -1 ../$ligname.$solvent.aveField.tmp >> ../$solvent.replicate.tmp # Output the value to calculate average of averages

                # Calculate the standard deviation
                awk '{x+=$4;y+=$4^2}END{print sqrt(y/NR-(x/NR)^2)}' counted_lines.txt > ../$ligname.$solvent.stdDev.tmp
                head -1 ../$ligname.$solvent.stdDev.tmp >> ../$solvent.stdDev.replicate.tmp # Output the value to calculate average of stDev

                # Put all columns together
                cd $workingDirectory/$runDirectory/$ligname/
                paste $ligname.$solvent.aveField.tmp $ligname.$solvent.stdDev.tmp $ligname.$solvent.noremovedlines.tmp > $solvreplica.tmp

                # Put all replicates of the same solvent together
                paste combinedsolvents.tmp $solvreplica.tmp > combination.tmp
		rm $solvreplica.tmp
                cp combination.tmp combinedsolvents.tmp

            done # End repeat for all replicates

            cp combinedsolvents.tmp $ligname.$solvent.tmp

            # Calculate average of replicates
            sed -i '1d' $solvent.replicate.tmp          # Remove the top blank line
            awk '{ total += $1 } END { print total/NR }' $solvent.replicate.tmp > $solvent.replicate.aveField.tmp

            sed -i '1d' $solvent.stdDev.replicate.tmp          # Remove the top blank line
            awk '{ total += $1 } END { print total/NR }' $solvent.stdDev.replicate.tmp > $solvent.replicate.aveStdDev.tmp

            # Calculate the standard error of the replicates
            awk '{x+=$1;y+=$1^2}END{print (sqrt(y/NR-(x/NR)^2))/sqrt(NR)}' $solvent.replicate.tmp > $solvent.replicate.stdErr.tmp

            # Count the number of data points
            awk '{SUM+=$1}END{print SUM}' $solvent.counted.tmp > $solvent.replicate.counted.tmp

            # Paste the average and standard error into separate file
            paste $ligname.$solvent.name.tmp $solvent.replicate.aveField.tmp $solvent.replicate.stdErr.tmp $solvent.replicate.aveStdDev.tmp $solvent.replicate.counted.tmp > $ligname.stdErr.tmp
            cat final.tmp $ligname.stdErr.tmp > combined.final.tmp
            cp combined.final.tmp final.tmp
            
            # Put all solvents together
            cat combinedfields.tmp $ligname.$solvent.tmp > combined.tmp
            cp combined.tmp combinedfields.tmp

        done # End repeat for all solvents

        cp combinedfields.tmp $ligname.$boxsize.fields.log
        cp final.tmp $ligname.$boxsize.aveFields.log
        cat $workingDirectory/summary.tmp $ligname.$boxsize.aveFields.log > $workingDirectory/combinedave.tmp
        cp $workingDirectory/combinedave.tmp $workingDirectory/summary.tmp

        # Tidy up 
        rm *.tmp

        cd $workingDirectory/$runDirectory/$ligname

    done # End repeat for all boxsizes

    cd $workingDirectory/$runDirectory

done
cd $workingDirectory
mv summary.tmp $identifier.log
rm *.tmp
