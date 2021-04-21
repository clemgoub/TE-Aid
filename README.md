# consensus2genome
consensus2genome (c2g) blasts any transposable element (TE) consensus sequence (nucleotides, fasta) to a reference genome and outputs a mapping graph displaying blastn hits divergence to their consensus as horizontal lines (y axis) along the consensus (x axis). It also plots the overall consensus coverage after piling-up all copies.

The function performs the **blastn** automatically through the system and print a customizable graph with **R**.

The consensus size, number of fragments (hits) and full length copies (according to user-defined threshold) found are automatically printed on the graph.

support: click the "issues" tab on github or [email me](mailto:goubert.clement@gmail.com)

<img src=https://github.com/clemgoub/consensus2genome/blob/master/Example/cons2gen.jpeg width="900">

*************
changelog v2
- add shell wrapper, coverage on second graph
changelog v1.1
- added the coverage curve (right y axix)
*************

## Dependencies
- [R (Rscript)](https://cran.r-project.org/mirrors.html)
- [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

consensus2genome calls blastn and R from the command line with 'blastn' and 'Rscript' commands. Both `blastn` and `Rscript` must be in the user path (usually the case following the default install of these progams). 
If not, you need to locate the executables location and add them to your local path before using c2g: 
```
export PATH="/path/to/blast/bins/folder/:$PATH"` 
export PATH="/path/to/R/bins/folder/:$PATH"` 
```

## Using the shell wrapper (recommended)

### Install from github

```
git clone https://github.com/clemgoub/consensus2genome.git
```

### Usage and option
#### Blastn databse
You will need to make a blastn database for your reference genome
```
makeblastdb -in genome.fa -out genome.fa -dbtype 'nucl'
```
### Usage
```
 cd your/path/to/consensus2genome/
 ./c2g.sh [-q|--query <query.TE.fa>] [-d|--blast-database <genome.fa>] [options]
```

#### Arguments

Mendatory arguments:
```
    -q, --query                   TE consensus to blast (fasta file)
    -d, --blast-database          Reference genome in blastn format database (makeblastdb -in genome.fa -out genome.fa -dbtype 'nucl')
```
- Optional arguments:
```
    -h, --help                    show this help message and exit

    -o, --output                  output folder

    -e, --e-value                 blastn: e-value threshold to keep hit (default: 10e-8)
    -f, --full-length-threshold   blastn: min. proportion (hit_size)/(consensus_size) to be considered "full length" (0-1; default: 0.9)

    -a, --alpha                   graphical: transparency value for blastn hit segments (0-1; default 0.3)
    -F, --full-length-alpha       graphical: transparency value for full-length blastn hits segments (0-1; default 1)
    -y, --auto-y                  graphical: manual override for y axis max value (default: TRUE; otherwise: -y NUM)
```

Tutorial coming soon!

## Usage as R function

Originaly, c2g was simply a R function. Below is the original tutorial using the R console

### 1. Load the function in R
In R, simply copy and paste the content of consensus2genome.R into a R console and press Enter. Alternatively, you can 'source' the consensus2genome.R file:
```Rscript
source("consensus2genome.R")
```
The function is now loaded for your current R session.


### 2. Run consensus2genome
```Rscript
consensus2genome(query, db, FL_thresh, alpha, full_alpha, auto_y)
```
#### mendatory arguments
- **query** the path of your query file in fasta

- **db** the path of your blast database (nucleotide format)

#### optional arguments

- **evalue** The evalue threshold to keep a blastn hit, default = 10e-8

- **FL_thresh** Full-lenght threshold, in % of the consensus sequence. Will display in red the genomic hits >= to this threshold. Default = 90%

- **alpha** 0 to 1. Transparency of the hits displayed on the graphs (0 = invidible, 1 = dense), default = 0.3

- **full_alpha** 0 to 1. Transparency of the full lenght hits (according the FL_thresh) displayed on the graph (red) (0 = invidible, 1 = dense), default = 1

- **auto_y** T or 0 to N. auto adjustment of the y-axis. If true (default) the y-axis (% divergence of the hit to the consensus) is adjusted relative to the data. Can be manually adjusted by changing it with any value > 0 (in % divergence). 

### R Example

- In this example we are going to map the Gypsy-2 and Jockey elements of ***Drosophila melanogaster*** over the reference genome. We assume that we start from the main folder of consensus2genome package. The consensus sequence is located in the 'Example' folder, and we will need to download the ***D. melanogaster*** reference genome and make a blast database out of it. Let's go!

#### 1. Download the ***D. melanogaster*** genome and create a blast db
```shell
cd Example
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.17.fasta.gz
gunzip dmel-all-chromosome-r6.17.fasta.gz
makeblastdb -in dmel-all-chromosome-r6.17.fasta -out dmel-all-chromosome-r6.17.fasta -dbtype 'nucl'
```

#### 2. Open R (type 'R' in the terminal then press Enter), load consensus2genome and run the function
```Rscript
source("../consensus2genone.R") # load the function
consensus2genome("Gypsy2_DM.fasta", "dmel-all-chromosome-r6.17.fasta")
```
<img src=https://github.com/clemgoub/consensus2genome/blob/master/Example/Gypsy_example.jpeg width="550">

- As you can see, the graph tells you that the consensus is 7221 bp long, has 417 fragments (hits) on the reference genome and only one fragment is superior of equal to 90% of the consensus sequence (and it is displayed in red).

- Now lets play with the different options of the function, that can be useful with a more fuzy graph at a first glance. We are going to map Jockey, a recent LINE element. According to the evolutionaty biology of this TE family, we expect to see a lot of recent copies (little divergence), many full length copies, as well as the characteristic pattern of 5' truncation of the LINE retroelements.

```Rscript
consensus2genome("Jockey_DM.fasta", "dmel-all-chromosome-r6.17.fasta")
```

<img src=https://github.com/clemgoub/consensus2genome/blob/master/Example/Jockey-1.jpeg width="550">

```Rscript
consensus2genome("Jockey_DM.fasta", "dmel-all-chromosome-r6.17.fasta", full_alpha=0.2, alpha=0.2)
```
<img src=https://github.com/clemgoub/consensus2genome/blob/master/Example/Jockey-2.jpeg width="1100">

- As expected we observe a lot of recent copies, as well as 12 full lenght elements. However they seem to mask the other copies. In the same fashin, there are more hits on the 3' end than on the 5' end due to the 5' truncation. We can play with the transparency of the graph to have a better view.
What we've done? We have a little increased the transparency of the black lines (alpha=0.3 by default to alpha=0.2) and drasticly increased the one of the red 'full-length' hits (alpha = 1 by default to alpha = 0.2).

- Now we still don't really see well what is going on because most of the copies are recent, but the graph shows the whole range of divergence, including two older hits, up to 25%. We can "zoom" in to get a better view with the **auto_y** parameter

```Rscript
consensus2genome("Jockey_DM.fasta", "dmel-all-chromosome-r6.17.fasta", full_alpha=0.2, alpha=0.5, auto_y=3.5)
```
<img src=https://github.com/clemgoub/consensus2genome/blob/master/Example/Jockey-3.jpeg width="550">

- Note that I also adjusted here again **alpha** to 0.5 and **full_alpha** to 0.2 for a better render.

#### 3. Exporting the graph

```Rscript
pdf("Gypsy2.pdf", width = 8, height= 10)
consensus2genome("Gypsy2_DM.fasta", "dmel-all-chromosome-r6.17.fasta")
dev.off()
```
This will create a file called "Gypsy2.pdf" in your working directory (***i.e.*** 'Example' here). You can change this in with any path in the pdf("yourpath/file.pdf", with = , height = ) function of the code. The width and height can be left empty and are relative values. If you use large numbers for pdf, the text might appear very small. Alternatively, you can use the png() function that works the same to output you graph in png. However you should then adjust the dimention on pixels.

