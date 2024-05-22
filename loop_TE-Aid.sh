#! /bin/bash

cons=$(cat $1) # consensus headers
library=$2
genomes=$3

# get script launch dir, from https://stackoverflow.com/a/246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

for i in $cons
do
	outname=$(echo "$i" | sed 's/#/\t/g' | cut -f 1)
	$DIR/extractfasta.sh H $i $library > TE.fasta
	QHEAD=$(grep '>' TE.fasta | sed 's/>//g')
	$DIR/TE-Aid -q TE.fasta -g $genomes
	mv TE.fasta.c2g.pdf ${outname}.pdf
	rm TE.fasta.fai
done

# cleaning
echo "cleaning up..."
rm blastn.txt TE.fasta

echo ""
echo ""
echo " -------------"
echo "( FINISHED!!! )"
echo " -------------"
