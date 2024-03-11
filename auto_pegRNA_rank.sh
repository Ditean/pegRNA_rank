#!/bin/bash

# while IFS= read -r amplicon; do python3 rank.py --amplicon "$amplicon" --insert atgatcctgacgacggagaccgccgtcgtcgacaagcc; done < sequences.txt


amplicon_file=""
insert_file=""

# Parse command-line for user files
while [[ "$#" -gt 0 ]]
do
	case "$1" in
		-a|--amplicon)
			amplicon_file="$2"
			shift 2
			;;

		-i|--insert)
			insert_file="$2"
			shift 2
			;;

		*)
			echo "Unknown option: $1"
			exit 1
			;;
	esac
done

# Check if required options are available
if [ -z "$amplicon_file" ] || [ -z "$insert_file" ]
then
	echo "Usage: $0 --amplicon <amplicon_file> --insert <insert_file>"
	exit 1
fi

# Check if amplicon file exists
if [ ! -f "$amplicon_file" ]
then
	echo "Error: Amplicon file '$amplicon_file' not found."
	exit 1
fi

# Check if insert file exists
if [ ! -f "$insert_file" ]
then
	echo "Error: Insert file '$insert_file' not found."
	exit 1
fi


# Loop through each line of amplicon file
while IFS= read -r amplicon
do

	# Loop through each line in insert file
	while IFS= read -r insert
	do
		# echo "python3 rank.py --amplicon $amplicon --insert 
$insert"
		python3 ~/Desktop/pegRNA_rank/rank.py --amplicon 
"$amplicon" --insert "$insert"

	done < "$insert_file"

done < "$amplicon_file"