#! /bin/bash
#######################################################################################
### extractfasta.sh - V1 - Clement Goubert (2020) - goubert.clement@gmail.com       ###
### ------------------------------------------------------------------------------- ###
### This script extract a fasta sequence based on header or list of headers         ###
### USAGE:  ./extractfasta.sh $1<H/L: H=header given; L=list given>                 ###
###        $2<Header [string] / List of headers [file]> $3<Consensus.file.fasta>    ###
### query: path to query (fasta file)                                               ###
### db: path to blast db (blast formated nucleotide database)                       ###
#######################################################################################

# USAGE:
# ./extractfasta.sh $1<H/L: H=header given; L=list given> $2<Header [string] / List of headers [file]> $3<Consensus.file.fasta>

if [ $1  == "H" ]
then
#echo "H!"
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(echo "$2") $3
else
  #echo "L!"
  perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $2 $3
fi