## Tutorial

In this example we are going to some transposable elements of ***Drosophila melanogaster***. The consensus sequences for this tutorial are located in the `Example/` folder, and you will need to download the ***D. melanogaster*** reference genome (dm6) and make a blast database out of it. Let's go!

#### 1. Download the ***D. melanogaster*** genome and create a blast db

```shell
curl -o Example/dm6.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip Example/dm6.fa.gz
makeblastdb -in Example/dm6.fa -out Example/dm6.fa -dbtype 'nucl'
```

Two *D. melanogaster* TE consensus are present in the folder `Examples`

#### 2. Analyze *D. melanogaster* consensus TEs

Let's start with Jockey, a recent LINE element in the *D. melanogaster* genome

```shell
./TE-Aid ./TE-Aid -q Example/Jockey_DM.fasta -d Example/dm6.fa -o ../dm6example -m 500
```
<img src=https://github.com/clemgoub/TE-Aid/blob/master/Example/Jockey_new.jpeg width="550">

- As expected we observe a lot of recent copies, as well as 12 full lenght elements. However they seem to mask the other copies. In the same fashin, there are more hits on the 3' end than on the 5' end due to the 5' truncation. We can play with the transparency of the graph to have a better view.
What we've done? We have a little increased the transparency of the black lines (alpha=0.3 by default to alpha=0.2) and drasticly increased the one of the red 'full-length' hits (alpha = 1 by default to alpha = 0.2).



