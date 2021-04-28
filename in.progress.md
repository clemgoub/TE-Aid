## Tutorial

In this example we are going to map the Gypsy-2 and Jockey elements of ***Drosophila melanogaster*** over the reference genome. We assume that we start from the main folder of consensus2genome package. The consensus sequence is located in the 'Example' folder, and we will need to download the ***D. melanogaster*** reference genome and make a blast database out of it. Let's go!

#### 1. Download the ***D. melanogaster*** genome and create a blast db

```shell
curl -o Example/dm6.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip Example/dm6.fa.gz
makeblastdb -in Example/dm6.fa -out Example/dm6.fa -dbtype 'nucl'
```

Two *D. melanogaster* TE consensus are present in the folder `Examples`

#### 2. Analyze *D. melanogaster* consensus TEs

Let's start with Gypsy XXX, a famous LTR element.

```shell
./c2g.sh -q Examples/.fasta -d dm6.fa -o dm6_TEaid
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

