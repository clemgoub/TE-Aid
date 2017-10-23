# consensus2genome
consensus2genome is a R function that blast any TE consensus sequence (nucleotides, fasta) to a reference reference genome and then output a mapping graph which displays the blast hits along the consensus (x axis) according to their divergence to it (y axis).
The function performs automatically the blast through the system and print a customizable graph.

## Dependencies
- [R](https://cran.r-project.org/mirrors.html)
- [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

consensus2genome calls blastn from the command line as 'blastn', thus blastn must be in your path. If it is not the case, for example if you have installed it in your home directory, you can modify the code with the full path of your blastn executable ([email me](mailto:goubert.clement@gmail.com) if necessary)

## Usage

### 1. Load the function in R
In R, simply copy and paste the content of consensus2genome.R into a R console and press Enter. Alternatively, you can 'source' the consensus2genome.R file:
```R
source("consensus2genome.R")
```
The function is now loaded for your current R session.


### 2. Run consensus2genome
```
consensus2genome(query, db, FL_thresh, alpha, full_alpha, auto_y)
```
#### mendatory arguments
- **query** the path of your query file in fasta

- **db** the path of your blast database (nucleotide format)

#### optional arguments
**FL_thresh** Full-lenght threshold, in % of the consensus sequence. Will display in red the genomic hits >= to this threshold. Default = 90%

- **alpha** 0 to 1. Transparency of the hits displayed on the graphs (0 = invidible, 1 = dense), default = 0.3

- **full_alpha** 0 to 1. Transparency of the full lenght hits (according the FL_thresh) displayed on the graph (red) (0 = invidible, 1 = dense), default = 1

- **auto_y** T or 0 to N. auto adjustment of the y-axis. If true (default) the y-axis (% divergence of the hit to the consensus) is adjusted relative to the data. Can be manually adjusted by changing it with any value > 0 (in % divergence). 
