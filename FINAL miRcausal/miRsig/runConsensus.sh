#!/bin/bash

#*******************************************************************************
#Author: Joseph J. Nalluri
#Script: To implement consensus-based network inference based on 6 algorithms - 1. Pearson 2.Spearman 3. GENIE3 4. CLR 5. MRNETB 6. Distance correlation
#Date: 11/14/2016
# Usage: Keep the directory locations intact. 
# Run on Linux: ./runConsensus.sh [input_file]
# E.g.: ./runConsensus.sh input_file.txt
# *****************************************************************************

#Removing previous zipped folder, if generated
if [ -f Individual/NI_results/network_inference_results.zip ]
 then rm Individual_NI_results/network_inference_results.zip
fi

if [ -f "consensus.txt" ]
 then rm "consensus.txt"
fi

printf "\n\n"
printf "****************************************************************\n"
printf "This script will generate network based on consensus methodology\n"
printf "****************************************************************\n"
printf "\n"


FILEARG="$1" # Assigning the input file argument to variable
FILENAME=${FILEARG##*/} # Removing the '/' and '..' from the argument
FILE=${FILENAME%.*} # File will contain the input file name

# Pearson Correlation
python correlations.py "$1" pearson &
# Spearman Correlation
python correlations.py "$1" spearman &
# Kendall Correlation
python correlations.py "$1" kendall &
#CLR, MRNETB, GENIE3 and Distance Correlation script
Rscript CLR.r "$1" 1>> Log/script_execution_log.txt &
Rscript MRNETB.r "$1" 1>> Log/script_execution_log.txt &
Rscript GENIE.r "$1" 1>> Log/script_execution_log.txt &
Rscript distcorrelation.r "$1" 1>> Log/script_execution_log.txt &

wait

printf "Done with algorithms\n"
printf "Finding files...\n"
files=( pearson_predictions.txt spearman_predictions.txt kendall_predictions.txt CLR_predictions.txt MRNETB_predictions.txt GENIE3_predictions.txt DC_predictions.txt )
all_files_present=true

for i in "${files[@]}"
do
  if ! [ -f $i ]
   then
   $all_files_present=false
   printf "ERROR: Could not locate $i file\n"
   break
  fi
done

if [ $all_files_present = true ]
 then
 printf "All files found\n"
 zip Individual_NI_results/network_inference_results.zip pearson_predictions.txt spearman_predictions.txt kendall_predictions.txt CLR_predictions.txt MRNETB_predictions.txt GENIE3_predictions.txt DC_predictions.txt
 
 if [ -f Individual_NI_results/network_inference_results.zip ]
  then
  printf "Zipped folder generated... \nPerforming consensus...\n"
  java -Xmx200g -jar AverageRank/ConsensusNet.jar --cut -1 --file1 Individual_NI_results/network_inference_results.zip >> Log/script_execution_log.txt
  
  if [ -f consensus.txt ]
   then 
   printf "\n***********************************************************************\n"
   printf "SUCCESS!! Final consensus-based network generated as consensus.txt file\n"
   
   else printf "\n\nERROR: consensus.txt not found\n\n\n"
  fi
  
  else printf "\n\nERROR: zipped folder not found\n\n\n"
 fi
 else printf "\n\nERROR: all files not found\n\n\n"
fi
 

if [ -f Individual_NI_results/network_inference_results.zip ]
 then printf "Zipped folder generated... \nPerforming consensus...\n"
 java -Xmx200g -jar AverageRank/ConsensusNet.jar --cut -1 --file1 Individual_NI_results/network_inference_results.zip >> Log/script_execution_log.txt
fi
  


#Removing individual algorithms result file since they are copied in network_inference_results.zip folder
rm pearson_predictions.txt
rm spearman_predictions.txt
rm kendall_predictions.txt
rm CLR_predictions.txt
rm MRNETB_predictions.txt
rm GENIE3_predictions.txt
rm DC_predictions.txt


printf "\n\n"
