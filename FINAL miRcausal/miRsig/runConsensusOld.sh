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

# Pearson Correlation
printf "1. Executing Pearson predictions...\n"
java -jar BasicCorrelation/Correlation.jar --data "$1" --cut -1 --cor pearson > Log/script_execution_log.txt

FILEARG="$1" # Assigning the input file argument to variable
FILENAME=${FILEARG##*/} # Removing the '/' and '..' from the argument
FILE=${FILENAME%.*} # File will contain the input file name

if [ -f "$FILE"_Correlation_predictions.txt ]
 then
  mv "$FILE"_Correlation_predictions.txt pearson_predictions.txt
  printf "...Pearson predictions output file generated.\n"
fi

# Spearman Correlation
printf "2. Executing Spearman predictions...\n"

java -jar BasicCorrelation/Correlation.jar --data "$1" --cut -1 --cor spearman 1>> Log/script_execution_log.txt

if [ -f "$FILE"_Correlation_predictions.txt ]
 then
  mv "$FILE"_Correlation_predictions.txt spearman_predictions.txt
  printf "...Spearman predictions output file generated.\n"
fi

#CLR, MRNETB, GENIE3 and Distance Correlation script
Rscript networkinference_mrnetB_clr_genie.r "$1" 1>> Log/script_execution_log.txt

printf "3. Executing CLR predictions...\n"
if [ -f CLR_predictions.txt ]
 then
 printf "... CLR predictions output file generated.\n"
 printf "4. Executing MRNETB predictions...\n"
  if [ -f MRNETB_predictions.txt ]
   then
   printf "... MRNETB predictions output file generated.\n"
   printf "5. Executing GENIE3 predictions...\n"
     if [ -f GENIE3_predictions.txt ]
       then
       printf "... GENIE3 predictions output file generated.\n"
       printf "6. Distance Correlation predictions...\n"
         if [ -f DC_predictions.txt ]
         then 
          printf "... Distance Correlation predictions output file generated\n"
          printf  "\nZipping all output files into network_inference_results folder... \n"
          zip Individual_NI_results/network_inference_results.zip pearson_predictions.txt spearman_predictions.txt CLR_predictions.txt MRNETB_predictions.txt GENIE3_predictions.txt DC_predictions.txt
          else printf "ERROR: Could not zip all the files...\n" 
         fi
       else printf "ERROR: Could not locate MRNETB_predictions.txt file\n"
     fi
   else printf "ERROR: Could not locate GENIE3_predictions.txt file\n"
  fi
 else printf "ERROR: Could not locate CLR_predictions.txt file\n"
fi

if [ -f Individual_NI_results/network_inference_results.zip ]
 then printf "Zipped folder generated... \nPerforming consensus...\n"
 java -Xmx200g -jar AverageRank/ConsensusNet.jar --cut -1 --file1 Individual_NI_results/network_inference_results.zip >> Log/script_execution_log.txt
fi
  
if [ -f consensus.txt ]
 then 
 printf "\n***********************************************************************\n"
 printf "SUCCESS!! Final consensus-based network generated as consensus.txt file\n"
fi

#Removing individual algorithms result file since they are copied in network_inference_results.zip folder
rm pearson_predictions.txt
#rm spearman_predictions.txt
rm CLR_predictions.txt
rm MRNETB_predictions.txt
rm GENIE3_predictions.txt
rm DC_predictions.txt


printf "\n\n"
