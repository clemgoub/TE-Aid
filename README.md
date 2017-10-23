# consensus2genome
consensus2genome is a R function that blast any TE consensus sequence (nucleotides, fasta) to a reference reference genome and output a mapping graph displaying the blast hits relative to the consensus.
The function performs automatically the blast through the system and print a customizable graph.

## dependancies
[R](https://cran.r-project.org/mirrors.html)
[blastn](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
consensus2genome calls blastn from the command line as 'blastn', thus blastn must be in your path. If it is not the case, for example if you don't have root access, you can modify the code with the full path of your blastn executable ([email me](goubert.clement@gmail.com) if necessary)

## usage
