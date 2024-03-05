#! /bin/Rscript
#######################################################################################
### consensus2genome - v2.1 - Artem Ilin (2023) - liartom2@gmail.com      ###
### ------------------------------------------------------------------------------- ###
### This R function blasts a TE consensus against a reference genome and then plots  ###
### the genomic fragments found relative to the consensus sequence                  ###
### see https://github.com/clemgoub/consensus2genome for the full documentation     ###
### to use, copy and paste the following code into a R console                      ###
### USAGE: consensus2genome(query, db , ...)                                        ###
### query: path to query (fasta file)                                               ###
### db: path to blast db (blast formated nucleotide database)                       ###
#######################################################################################
# Changelog V2 --> V3 | 04.19.2021
# - Fit for shell wrapper, remove regions highlight from usused trim script
# Changelog V1 --> V2 | 03.13.2020
# - Add a second graph with suggested cut points
# 
# V2=alpha | not for release
# TO FIX: trace the coverage on the left graph
# TO DO: convert breakpoints to bed for getfasta



consensus2genome=function(query=NULL, db=NULL, evalue=10e-8,
                          FL_thresh=0.9, alpha=0.3, full_alpha=1,
                          auto_y=T, wdir=NULL){
  if(!requireNamespace("Biostrings", quietly = T)){
    stop("Biostrings not installed")
  } else if (!requireNamespace("Rcpp", quietly = T)){
    stop("Rcpp not installed, quitting...")
  } else {
  if(is.null(query)){print('query not specified')}
  if(is.null(db)){print('db not specified')}
  #perform the blast
  getwd()
  blast=read.table(text=system(paste("blastn -query", query, "-db", db , "-evalue", evalue, "-outfmt 6 | sed 's/#/-/g'"), intern = TRUE), colClasses = c("character", "character", "numeric", rep("integer", 7), "numeric", "numeric"))
  colnames(blast) <- c("query", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "score")
  # Load reduce_blast function
  Rcpp::sourceCpp(file.path(wdir, "reduce.cpp"))
  blast <- reduce_blastn_results(blast)
  write.table(blast, file = "blastn.txt",  quote = F, row.names = F, col.names = F)
  #TE consensus size
  cons_len <- Biostrings::fasta.seqlengths(query)
  print(paste("consensus length:", cons_len, "bp", sep = " "))
  #list of almost full length fragments
  full=blast[abs(blast$qend-blast$qstart) >= FL_thresh*as.numeric(cons_len),]
  #graph
  if(auto_y == T){
    #par(mar=c(2,2,2,2))
    plot(range(0, cons_len), range(0, max(100-blast$pident)), type = "n", main=paste("TE: ", as.character(blast[1,1]), "\n consensus size: ", as.character(cons_len), "bp; fragments: ", as.character(length(blast$query)), "; full length: ", as.character(length(full$query))," (>=",as.character(as.numeric(FL_thresh)*cons_len),"bp)", sep = ""), cex.main = 0.9, xlab = "TE consensus (bp)", ylab = "divergence to consensus (%)")
    for(i in 1:nrow(blast)){
      segments(blast$qstart[i], 100-blast$pident[i], blast$qend[i], 100-blast$pident[i], col=rgb(0,0,0,alpha=alpha))
    }
    for(i in 1:length(blast$query)){
      segments(full$qstart[i], 100-full$pident[i], full$qend[i], 100-full$pident[i], col=rgb(1,0,0, alpha=full_alpha), lwd=1.5)
    }}
  else{
    #par(mar=c(2,2,2,2))
    plot(range(0, cons_len), range(0, auto_y), type = "n", main=paste("TE: ", as.character(blast[1,1]), "\n consensus size: ", as.character(cons_len), "bp; fragments: ", as.character(length(blast$query)), "; full length: ", as.character(length(full$query))," (>=",as.character(as.numeric(FL_thresh)*cons_len),"bp)", sep = ""), cex.main = 0.9, xlab = "TE consensus (bp)", ylab = "divergence to consensus (%)")
    for(i in 1:length(blast$query)){
      segments(blast$qstart[i], 100-blast$pident[i], blast$qend[i], 100-blast$pident[i], col=rgb(0,0,0,alpha=alpha))
    }
    for(i in 1:length(blast$query)){
      segments(full$qstart[i], 100-full$pident[i], full$qend[i], 100-full$pident[i], col=rgb(1,0,0, alpha=full_alpha), lwd=1.5)
    }}
  
#make the coverage matrix and graph
coverage=matrix(rep(0, (length(blast$query)*as.numeric(cons_len))), byrow = T, ncol = as.numeric(cons_len))
# print(dim(coverage))
for(i in 1:nrow(blast)){
    coverage[i, blast$qstart[i]:blast$qend[i]] <- 1
}

    # TO FIX: trace the coverage on the left graph
    #points(colSums(coverage), type='l', axes = F, ylab = NA, xlab = NA, col=covcol, ylim = c(0, max(colSums(coverage))))
    #axis(side = 4)
    #mtext(side = 4, line = 3, 'consensus coverage (bp)')
    
    ## import removator
    
    removator<-function(covM){
      covMT <- as.data.frame(covM)
      covMT$bp <- rownames(covMT)
      plot(covM, type = "l", main = "TE consensus genomic coverage", xlab = "TE consensus (bp)", ylab = "coverage (bp)")
     } # removator function
    
    removator(colSums(coverage)) # makes the second graph
  }
}
