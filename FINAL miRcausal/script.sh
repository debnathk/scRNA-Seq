
# INPUT: input file, score, reset value, percent to keep
rm -r intermediate_outputs
rm ../miRsig_pipeline/Consensus_Alcohol_Dataset/Consensus_algos/consensus.txt
rm -r other_files

mkdir intermediate_outputs
mkdir other_files
FILE="$1"
what_to_run="$2"
score="$3"
reset_val="$4"
percent="$5"
if [ -f "influential_mirna.txt" ]; then rm influential_mirna.txt; fi
if [ -f "mirsig_network.txt" ] && [ "$what_to_run" -ne "2" ]; then rm mirsig_network.txt; fi
echo "File: " $FILE
echo "Cutoff score: " $score
echo "Reset value: " $reset_val
echo "Decimal percent to keep: " $percent
#:'
case $what_to_run in 
	"1")
		echo "Running only miRsig"
		;;
	"2")
		echo "Running only miRfluence"
		;;
	"3")
		echo "Running miRsig and miRfluence"
		;;
esac


if [ -f "$FILE" ]
 then
	if [ "$what_to_run" == "1" ] || [ "$what_to_run" == "3" ]
	 then
		echo "Starting miRsig processing"
		set -e
		python preprocess.py "$FILE"
		
		cd miRsig/

		./runConsensus.sh  ../other_files/cleaned_matrix.txt
		cd ../
		echo "Done with miRsig"
		if [ "$what_to_run" == "1" ]
		 then 
			cp miRsig/consensus.txt mirsig_network.txt
			python map_back_mirsig.py
		 else
			cp miRsig/consensus.txt intermediate_outputs/mirsig_network.txt
		fi
		echo "Copied consensus.txt"
	fi
	if [ "$what_to_run" == "2" ] || [ "$what_to_run" == "3" ]
	 then
		echo "Starting miRfluence processing"
		if [ "$what_to_run" == "2" ]
		 then
			python process_for_mirfluence.py "$FILE" "$score" "$reset_val"
		 else
			python process_for_mirfluence.py intermediate_outputs/mirsig_network.txt "$score" "$reset_val"
		fi
		
		echo "Processed for miRfluence"
		if [ -f "miRfluence/coverage_inference.txt" ]
			then
				rm miRfluence/coverage_inference.txt
				echo "Removed old miRfluence output"
		fi
		
		if [ "$what_to_run" == 3 ]
		  then
			first_line="probGraphFile : ../other_files/filtered_consensus.txt"
		  else
		    first_line="probGraphFile : $FILE"
		fi
		cd miRfluence

		sed -i "s|probGraphFile : .*|$first_line|" input.txt
		sed -i "s|outdir : .*|outdir : ../intermediate_outputs/|" input.txt

		
		echo "Set up miRfluence"

		echo "Running miRfluence..."
		./InfluenceModels -c input.txt
		cd ../
		if [ -f "miRfluence/coverage_inference.txt" ]
			then
				echo "Ran miRfluence"
				cp miRfluence/coverage_inference.txt intermediate_outputs/mirfluence_output.txt
				echo "Copied mirfluence"
				python filter_mirfluence_output.py intermediate_outputs/mirfluence_output.txt "$percent" 
				echo;echo
				echo "Top few influential miRNA:"
				head -10 influential_mirna.txt
				echo
				echo "**Full list in influential_mirna.txt**"
				echo;echo
			else
				echo "miRfluence Failed"
		fi
	fi
 else echo "File not found."
fi
#'