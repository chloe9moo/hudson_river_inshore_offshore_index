#!/bin/bash
##script to copy and rerun vast models for multiple species/surveys##

MAXJOBS=4

#file variables
RUNDIR="~/Documents/projects/hudson_river"
ORIGFILE="11_vast_cov_runs.R"
RUNFILENAME="cov_comp"
RUNDATE=$(date +%Y%m%d)

#in script variables
#SPECIES2RUN="Striped Bass"
#SURVEY2RUN="all"
SPECIES2RUN=("Alewife" "Striped Bass" "American Shad" "Blueback Herring")
SURVEY2RUN=("all" "fjs" "bss")
NKNOTS=200 #(50 100 150 200 250 300 500 800 1000)
BIASCORRECT="TRUE"
CALCRANGE="FALSE"
CALCEFFAREA="FALSE"
SKIPANYCATCHRUNS="FALSE"
SKIPCATCHVALS=""
SKIPANYDENSRUNS="FALSE"
SKIPDENSVALS=""

echo "Starting multi model run for ${RUNFILENAME}"
echo "Species included are ${SPECIES2RUN[*]}"
echo "Surveys included are ${SURVEY2RUN[*]}"
echo "Knots included are ${NKNOTS[*]}"
echo " "

#prep code for run
mkdir "wd_${RUNFILENAME}_${RUNDATE}"
cp $ORIGFILE "wd_${RUNFILENAME}_${RUNDATE}"
#head -49 $ORIGFILE >> "wd_${RUNFILENAME}_${RUNDATE}/${ORIGFILE}" ##testing
cd "wd_${RUNFILENAME}_${RUNDATE}" 

for spp in "${SPECIES2RUN[@]}"; do
  for surv in "${SURVEY2RUN[@]}"; do
    for k in "${NKNOTS[@]}"; do
    
    	SAVENAME="${RUNFILENAME}_${spp// /_}_${surv}_${k}"
   	cp $ORIGFILE ${SAVENAME}.R
   	
  	#set working directory for R script
	sed -i "s#PATH <- getwd()#PATH <- \"${RUNDIR}\"#" "${SAVENAME}.R"

	#replace species to run in R script
	sed -i 's/\(spp = "\)[^"]*\(".*\)/\1'"${spp}"'\2/' ${SAVENAME}.R	
	
	#replace survey to run in R script
	sed -i 's/\(survey2run = "\)[^"]*\(".*\)/\1'"${surv}"'\2/' ${SAVENAME}.R
	
	#replace knots num. to run in R script
	sed -i "s/^\(num_knots *= *\).*/\1${k}/" ${SAVENAME}.R
	
	#set in script variables
	#sed -i "s/^\(num_knots *= *\).*/\1${NKNOTS}/" ${SAVENAME}.R
	sed -i "s/^\(do.bias.correct *= *\).*/\1${BIASCORRECT}/" ${SAVENAME}.R
	sed -i "s/^\(calc_range_in_models *= *\).*/\1${CALCRANGE}/" ${SAVENAME}.R
	sed -i "s/^\(calc_eff_area_in_models *= *\).*/\1${CALCEFFAREA}/" ${SAVENAME}.R
	
	#set whether to skip any of the covariate combos
	sed -i "s/skip.catch = .*/skip.catch = ${SKIPANYCATCHRUNS}/" ${SAVENAME}.R
	if [ "$SKIPANYCATCHRUNS" == "TRUE" ]; then 
		echo "Skipping catch runs: ${SKIPCATCHVALS}"
	  	sed -E -i "s|i %in% c\\([^)]*\\)|i %in% c(${SKIPCATCHVALS})|" "${SAVENAME}.R"
  	fi
  
  	sed -i "s/skip.dens = .*/skip.dens = ${SKIPANYDENSRUNS}/" ${SAVENAME}.R
	if [ "$SKIPANYDENSRUNS" == "TRUE" ]; then
	  	echo "Skipping dens runs: ${SKIPDENSVALS}"
	  	sed -E -i "s|y %in% c\\([^)]*\\)|y %in% c(${SKIPDENSVALS})|" "${SAVENAME}.R"
  	fi


	#set up output log
	{
	     	echo "Run started at: $(date)"
		echo "Species run =                   ${spp}"
  		echo "Survey run  =                   ${surv}"
  		echo "Knots used  =                   ${k}"
		echo "Bias corrected?                 ${BIASCORRECT}"
  		echo "Range calculated?               ${CALCRANGE}"
  		echo "Effective area occ. calculated? ${CALCEFFAREA}"
  		echo "----------------------------------------------"
  		echo "----------------------------------------------"
		echo " "
		echo " "
	} > "out_${SAVENAME}.txt"

	
	#run R script for species x survey combo
	(
	 echo "Starting run for ${SAVENAME} at $(date +%H:%M)"
	 Rscript --vanilla ./${SAVENAME}.R >> "out_${SAVENAME}.txt" 2>&1
	 echo "Run completed for ${SAVENAME} at $(date +%H:%M)"
	 echo " "
 	) &

	#limit jobs
	while (( $(jobs -r | wc -l) >= MAXJOBS )); do
		sleep 5
	done
	
		done
	done
done

#make sure all finish
wait
echo "All runs completed at $(date +%H:%M:%S)"
