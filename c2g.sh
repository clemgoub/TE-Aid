#! /bin/bash

################################################
# c2g - consensus to genome R function wrapper #
################################################
#
# V.0 | 04.19.21 - first version of the wrapper for test 
#
# Author: Clément Goubert - goubert.clement@gmail.com

#########################################################################################
#### PARSER: from https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f ####
##
#

# usage function
function usage()
{
   cat << HEREDOC

   Usage: c2g [-q|--query FILE] [-d|--blast-database FILE] [options]

   mendatory arguments:
    
    -q, --query                   TE consensus to blast (fasta file)
    -d, --blast-database          Reference genome in blastn format database (makeblastdb -in genome.fa -out genome.fa -dbtype 'nucl')

   optional arguments:
    
    -h, --help                    show this help message and exit
    
    -o, --output                  output folder
    
    -e, --e-value                 blastn: e-value threshold to keep hit (default: 10e-8)
    -f, --full-length-threshold   blastn: min. proportion (hit_size)/(consensus_size) to be considered "full length" (0-1; default: 0.9)
    
    -m, --min-orf                 getorf: minimum ORF size (in bp)

    -a, --alpha                   graphical: transparency value for blastn hit (0-1; default 0.3)
    -F, --full-length-alpha       graphical: transparency value for full-length blastn hits (0-1; default 1)
    -y, --auto-y                  graphical: manual override for y lims (default: TRUE; otherwise: -y NUM)

HEREDOC
} 

# if no parameter given, output help and qui
if [[ $# -eq 0 ]] ; then
    echo '   **********************************'
    echo '   Error! No mendatory argument given'
    echo '   **********************************'
    usage
    exit 0
fi

# parameters parser
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -q|--query)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        QUERY=$2
        shift 2
      else
        echo "Error: no query (-q) provided (fasta)" >&2
        usage
        exit 1
      fi
      ;;
	-d|--blast-database)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        GENOME_DB=$2
        shift 2
      else
        echo "Error: No genome datablase provided (blastn database)" >&2
        usage
        exit 1
      fi
      ;;    
 	-e|--e-value)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        EVALUE=$2
        shift 2
      fi
      ;;
  	-f|--full-length-threshold)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        FL=$2
        shift 2
      fi
      ;;
    -a|--alpha)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        ALPHA=$2
        shift 2
      fi
      ;;
    -F | --full-length-alpha)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        FULL_ALPHA=$2
        shift 2
      fi
      ;;   
     -y | --auto-y)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        AUTO_Y=$2
        shift 2
      fi
      ;;  
      -m | --min-orf)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        MINORF=$2
        shift 2
      fi
      ;; 
    -h | --help)
	  usage
	  exit 1
	  ;;
	-o | --output)
	  if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
	    OUTPUT=$2
	    shift 2  
	  fi
      ;;   	
    -*|--*=) # unsupported flags
      echo "Error: Unsupported argument $1" >&2
      usage
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac # <- end of case
done
# set positional arguments in their proper place
eval set -- "$PARAMS"


#########################################################################################
#### MAIN:
##
#

# asign default value and print parameters
TENAME="$(echo "$QUERY" | sed 's/\//\t/g' | awk '{print $NF}')"
EVALUE="${EVALUE:-10e-8}"
FL="${FL:-0.9}"
ALPHA="${ALPHA:-0.3}"
FULL_ALPHA="${FULL_ALPHA:-0.9}"
AUTO_Y="${AUTO_Y:-TRUE}"
OUTPUT="${OUTPUT:-.}"
MINORF="${MINORF:-300}"

# param check
echo "query:                  $QUERY"
echo "genome db:              $GENOME_DB"
echo "e-value:                $EVALUE"
echo "full length min ratio:  $FL"
echo "hits transparency:      $ALPHA"
echo "full length hits trsp.: $FULL_ALPHA"


## run script
# create output dir if non-existent
mkdir -p $OUTPUT
# create dotmatcher graph in outputs
dotmatcher -asequence $QUERY \
           -bsequence $QUERY \
           -graph png \
           -windowsize 25 \
           -threshold 50 \
           -goutfile $OUTPUT/$TENAME"".dotmatcher.png
# make temporary database for self-dotplot
makeblastdb -in $QUERY -out $OUTPUT/TE.db -dbtype 'nucl'
# run getorf
getorf -sequence $QUERY --outseq $OUTPUT/TE.orfs -minsize $MINORF
grep '>' $OUTPUT/TE.orfs | awk '{print $2"\t"$4}' | sed 's/\[//g;s/\]//g' > $OUTPUT/TE.orfs.R
# run R script with user-defined parameters
Rscript Run-c2g.R $QUERY $GENOME_DB $EVALUE $FL $ALPHA $FULL_ALPHA $AUTO_Y $OUTPUT $OUTPUT/TE.db $OUTPUT/TE.orfs.R $MINORF
# clean-up
rm $OUTPUT/TE.db*

echo "Done! The graph (.pdf) can be found in the output folder: $OUTPUT"