#!/bin/bash
##script to read in vast models within a single species x survey and compare diagnostics##

#file variables
RUNDIR="~/Documents/projects/hudson_river"
ORIGFILE="14_vast_model_comp.R"
RUNFILENAME="vast_model_comp"
RUNDATE=$(date +%Y%m%d)

#in script variables
#SPECIES2RUN="Striped Bass"
#SURVEY2RUN="all"
SPECIES2RUN=("Alewife") #"Striped Bass") #"American Eel" "Spottail Shiner")
SURVEY2RUN=("all" "fjs" "bss")
DATE2CHECK=("20250530" "20250531" "20250601")

#prep code for run
mkdir "wd_${RUNFILENAME}_${RUNDATE}"
cp $ORIGFILE "wd_${RUNFILENAME}_${RUNDATE}"
#head -49 $ORIGFILE >> "wd_${RUNFILENAME}_${RUNDATE}/${ORIGFILE}" ##testing
cd "wd_${RUNFILENAME}_${RUNDATE}" 

for spp in "${SPECIES2RUN[@]}"; do
	for surv in "${SURVEY2RUN[@]}"; do
		for datecheck in "${DATE2CHECK[@]}"; do
			
			SAVENAME="${RUNFILENAME}_${spp// /_}_${surv}_${datecheck}"
   			cp $ORIGFILE ${SAVENAME}.R
	
		  	#set working directory for R script
	  		sed -i "s#PATH <- getwd()#PATH <- \"${RUNDIR}\"#" "${SAVENAME}.R"

	  		#replace species to run in R script
	  		sed -i 's/\(spp = "\)[^"]*\(".*\)/\1'"${spp}"'\2/' ${SAVENAME}.R	
	
	  		#replace survey to run in R script
	  		sed -i 's/\(survey2check = "\)[^"]*\(".*\)/\1'"${surv}"'\2/' ${SAVENAME}.R

	  		#replace file save date string in R script
	  		sed -i 's/\(date2check = "\)[^"]*\(".*\)/\1'"${datecheck}"'\2/' ${SAVENAME}.R

  			#set up output log
  			{
  	  			echo "Run started at: $(date)"
  				echo "Species run =  ${spp}"
    				echo "Survey run  =  ${surv}"
   				echo "Date checked = ${datecheck}" 
      				echo "----------------------------------------------"
    				echo "----------------------------------------------"
  				echo " "
  				echo " "
  			} > "out_${SAVENAME}.txt"

			#run R script for species x survey combo
			echo "Starting run for ${SAVENAME} at $(date +%H:%M)"
			Rscript --vanilla ./${SAVENAME}.R >> "out_${SAVENAME}.txt" 2>&1
			echo "Run completed for ${SAVENAME} at $(date +%H:%M)"
			echo " "
			wait

	
		done
	done
done


