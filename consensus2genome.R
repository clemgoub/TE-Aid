#! /bin/Rscript
#######################################################################################
### consensus2genome - v2 - Clement Goubert (2020) - goubert.clement@gmail.com      ###
### consensus2genome - v2.1 - Artem Ilin (2023) - liartom2@gmail.com                ###
### ------------------------------------------------------------------------------- ###
### This R function blasts a TE consensus against a reference genome and then plots ###
### the genomic fragments found relative to the consensus sequence                  ###
### see https://github.com/clemgoub/consensus2genome for the full documentation     ###
### to use, copy and paste the following code into a R console                      ###
### USAGE: consensus2genome(query, db , ...)                                        ###
### query: path to query (fasta file)                                               ###
### db: path to blast db (blast formated nucleotide database)                       ###
#######################################################################################
# Changelog V2 --> V3 | 04.19.2021
# - Fit for shell wrapper, remove regions highlight from used trim script
# Changelog V1 --> V2 | 03.13.2020
# - Add a second graph with suggested cut points
# 
# V2=alpha | not for release
# TO FIX: trace the coverage on the left graph
# TO DO: convert breakpoints to bed for getfasta



consensus2genome=function(query=NULL, db=NULL, evalue=10e-8,
                          FL_thresh=0.9, alpha=0.3, full_alpha=1,
                          auto_y=T, wdir=NULL, nored=NULL){
  if(!requireNamespace("Biostrings", quietly = T)){
    stop("Biostrings not installed")
  } else {
  if(is.null(query)){print('query not specified')}
  if(is.null(db)){print('db not specified')}
  #perform the blast
  getwd()
  blast=read.table(text=system(paste("blastn -query", query, "-db", db , "-evalue", evalue, "-outfmt 6 | sed 's/#/-/g'"), intern = TRUE), colClasses = c("character", "character", "numeric", rep("integer", 7), "numeric", "numeric"))
  colnames(blast) <- c("query", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "score")
  #TE consensus size
  cons_len <- Biostrings::fasta.seqlengths(query)
  print(paste("consensus length:", cons_len, "bp", sep = " "))
  
  if (nored == "TRUE"){
    if (!requireNamespace("Rcpp", quietly = T)){
      stop("Rcpp not installed, quitting...")
    }
    # Load reduce_blast function
    Rcpp::sourceCpp(file.path(wdir, "reduce.cpp"))
    blast_red <- reduce_blastn_results(blast)
    # Export the reduced table
    write.table(blast_red, file = "blastn.txt",  quote = F, row.names = F, col.names = F)
    #list of almost full length fragments
    full=blast_red[abs(blast_red$qend-blast_red$qstart) >= FL_thresh*as.numeric(cons_len),]
    plot_main <- c(
      TE = as.character(blast_red[1,1]),
      ConsLength = as.character(cons_len),
      Fragments = length(blast_red$query),
      FullLength = length(full$query),
      LenThreshold = as.character(as.numeric(FL_thresh)*cons_len)
    )
  } else {
    write.table(blast, file = "blastn.txt",  quote = F, row.names = F, col.names = F)
    #list of almost full length fragments
    full=blast[abs(blast$qend-blast$qstart) >= FL_thresh*as.numeric(cons_len),]
    plot_main <- c(
      TE = as.character(blast[1,1]),
      ConsLength = as.character(cons_len),
      Fragments = length(blast$query),
      FullLength = length(full$query),
      LenThreshold = as.character(as.numeric(FL_thresh)*cons_len)
    )
  }
  
  
  #graph
  if(auto_y == T){
    #par(mar=c(2,2,2,2))
    plot(range(0, cons_len), range(0, max(100-blast$pident)), type = "n", main=paste("TE: ", plot_main[["TE"]], "\n consensus size: ", plot_main[["ConsLength"]], "bp; fragments: ", plot_main[["Fragments"]], "; full length: ", plot_main[["FullLength"]]," ( >= ", plot_main[["LenThreshold"]]," bp )", sep = ""), cex.main = 0.9, xlab = "TE consensus (bp)", ylab = "divergence to consensus (%)")
    for(i in 1:nrow(blast)){
      segments(blast$qstart[i], 100-blast$pident[i], blast$qend[i], 100-blast$pident[i], col=rgb(0,0,0,alpha=alpha))
    }
    for(i in 1:length(blast$query)){
      segments(full$qstart[i], 100-full$pident[i], full$qend[i], 100-full$pident[i], col=rgb(1,0,0, alpha=full_alpha), lwd=1.5)
    }}
  else{
    #par(mar=c(2,2,2,2))
    plot(range(0, cons_len), range(0, auto_y), type = "n", main=paste("TE: ", plot_main[["TE"]], "\n consensus size: ", plot_main[["ConsLength"]], "bp; fragments: ", plot_main[["Fragments"]], "; full length: ", plot_main[["FullLength"]]," ( >= ", plot_main[["LenThreshold"]]," bp )", sep = ""), cex.main = 0.9, xlab = "TE consensus (bp)", ylab = "divergence to consensus (%)")
    for(i in 1:length(blast$query)){
      segments(blast$qstart[i], 100-blast$pident[i], blast$qend[i], 100-blast$pident[i], col=rgb(0,0,0,alpha=alpha))
    }
    for(i in 1:length(blast$query)){
      segments(full$qstart[i], 100-full$pident[i], full$qend[i], 100-full$pident[i], col=rgb(1,0,0, alpha=full_alpha), lwd=1.5)
    }}
  
#make the coverage matrix and graph

coverage=matrix(rep(FALSE, length(blast$V1)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))
for(i in 1:length(blast$V1)){
    coverage[i,]<-c(rep(FALSE,blast$V7[i]-1),rep(TRUE,abs(blast$V8[i]-blast$V7[i])+1), rep(FALSE,as.numeric(cons_len)-blast$V8[i]))
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
