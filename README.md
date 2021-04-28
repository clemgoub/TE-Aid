# TE-Aid

TE-Aid is a program aimed to help the manual curation of transposable elements (TE). It inputs a TE consensus sequence (fasta format) and require a reference genome formated as a blastn database. Using `R` and `ncbi blast+`, TE-Aid produces 4 figures reporting (1) the genomic hits with divergence to consensus, (2) the genomic coverage of the consensus, (3) a self dot-plot and (4) a structure analysis including: TIR and LTR suggestions, open reading frames (ORFs) and TE protein hits annotation.

*include figure here* <img src=https://github.com/clemgoub/TE-Aid/blob/master/Example/TE1.jpeg width="900">

**Pipeline:**

- The TE (idealy, candidate consensus sequence) is searched against the provided reference genome with `blastn` 
	- Fig 1: genomic hits (horizontal lines) are represented relative to the query (TE consensus), the y axis represent the `blastn` divergence
	- Fig 2: pileup of the genomic hits relative to position along the query (TE consensus)
- The query is then blasted agaisnt itself in order to detect micro repeats and inversions (putative TIRs, LTRs)
	- Fig 3: self dot-plot and Fig 4 (top): TIR and LTR are suggested suggestions (colored arrows)
	- Bonus: a self dot-plot with `emboss dotmatcher` is also produced as an extra figure
- Putative ORFs are searched with `emboss getorf` and the peptides queried against a freely available TE protein databse (ref details to add)
	- Fig 4: ORFs (black rectangles: + orientation; red rectangles: - orientation), TE protein hits 

The consensus size, number of fragments (hits) and full length copies (according to user-defined threshold) are automatically printed on the graph.
If any ORFs and protein hits are found, their locations relative to the consensus are printed in the `stdout`

support: click the "issues" tab on github or [email me](mailto:goubert.clement@gmail.com)

TE-Aid comes from `consensus2genome` that will be soon deprecated

## Install

### Dependencies

- [R (Rscript)](https://cran.r-project.org/mirrors.html)
- [NCBI Blast suite](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

TE-Aid calls blastn and R from the command line with 'blastn', 'blastp', 'makeblastdb' and 'Rscript' commands. Both `blastn` and `Rscript` must be in the user path (usually the case following the default install of these progams). 
If not, you need to locate the executables location and add them to your local path before using c2g: 
```
export PATH="/path/to/blast/bins/folder/:$PATH"` 
export PATH="/path/to/R/bins/folder/:$PATH"` 
```

### Instal from github
```
git clone https://github.com/clemgoub/TE-Aid.git
```

## Usage and option

### Blastn databse

You will need to create a blastn database for your reference genome

```
makeblastdb -in genome.fa -out genome.fa -dbtype 'nucl'
```

### Usage

```
<user-path>/TE-Aid [-q|--query <query.TE.fa>] [-d|--blast-database <genome.fa>] [options]
```
replace `<user-path` with the path of the downloaded `TE-Aid` folder.

### Arguments

Mendatory arguments:
```
    -q, --query                   TE consensus to blast (fasta file)
    -d, --blast-database          Reference genome in blastn format database (makeblastdb -in genome.fa -out genome.fa -dbtype 'nucl')
```
- Optional arguments:
```
    -h, --help                    show this help message and exit

    -o, --output                  output folder. Default = "."

    -e, --e-value                 genome blastn: e-value threshold to keep hit (default: 10e-8)
    -f, --full-length-threshold   genome blastn: min. proportion (hit_size)/(consensus_size) to be considered "full length" (0-1; default: 0.9)

    -a, --alpha                   graphical: transparency value for blastn hit segments (0-1; default 0.3)
    -F, --full-length-alpha       graphical: transparency value for full-length blastn hits segments (0-1; default 1)
    -y, --auto-y                  graphical: manual override for y axis max value (default: TRUE; otherwise: -y NUM)
```

Tutorial coming soon!