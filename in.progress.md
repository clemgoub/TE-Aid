## Tutorial

In this example we are going to some transposable elements of *Drosophila melanogaster*. The consensus sequences for this tutorial are located in the `Example/` folder, and you will need to download the *D. melanogaster* reference genome (dm6) and make a blast database out of it. Let's go!

#### 1. Download the *D. melanogaster* genome

```shell
curl -o Example/dm6.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip Example/dm6.fa.gz
```
*D. melanogaster* TE consensus are present in the folder `Examples`

#### 2. Analyze the TE consensus

Let's start with Jockey, a recent **LINE** element in the *D. melanogaster* genome

```shell
./TE-Aid -q Example/Jockey_DM.fasta -g Example/dm6.fa -o ../dm6example
```
<img src=https://github.com/clemgoub/TE-Aid/blob/master/Example/Jockey.TEaid.png width="1024">

Next is Gypsy-2, from the **LTR** lineage

```shell
./TE-Aid -q Example/Gypsy2_DM.fasta -g Example/dm6.fa -o ../dm6example
```
<img src=https://github.com/clemgoub/TE-Aid/blob/master/Example/Gypsy.TEaid.png width="1024">
