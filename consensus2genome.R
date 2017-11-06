#! /bin/Rscript
#######################################################################################
### consensus2genome - v1.1 - Clement Goubert (2017) - goubert.clement@gmail.com    ###
### ------------------------------------------------------------------------------- ###
### This R function blast a TE consensus against a reference genome and then plots  ###
### the genomic fragments found relative to the consensus sequence                  ###
### see https://github.com/clemgoub/consensus2genome for the full documentation     ###
### to use, copy and paste the following code into a R console                      ###
### USAGE: consensus2genome(query, db , ...)                                        ###
### query: path to query (fasta file)                                               ###
### db: path to blast db (blast formated nucleotide database)                       ###
#######################################################################################

consensus2genome=function(RMfile=NULL, TE=NULL, FL_thresh=0.9, alpha=0.3, full_alpha=1, auto_y=T, covcol="blue"){
  if(is.null(RMfile)){print('RepeatMasker file not specified')}
  if(is.null(TE)){print('TE consensus to report not specified')}
  #perform the blast
  RM_cons=read.table(text=system(paste("awk ' NR > 3 {if ($9 != \"C\") {print $1\"\t\"$2\"\t\"$10\"\t\"$11\"\t\"$12\"\t\"$13\"\t\"$14} else {print $1\"\t\"$2\"\t\"$10\"\t\"$11\"\t\"$14\"\t\"$13\"\t\"$12}}'", RMfile, " | grep -w", TE, "| sed 's/(//g;s/)//g'"), intern = TRUE))
  #TE consensus size
  RMcons_len=RM_cons$V6[1]+RM_cons$V7[1]
  #list of almost full length fragments
  full=RM_cons[abs(RM_cons$V5-RM_cons$V6) >= FL_thresh*as.numeric(RMcons_len),]
  #graph
  if(auto_y == T){
    par(mar=c(5,5,5,5))
    plot(range(0, RMcons_len), range(0, max(RM_cons$V2)), type = "n", main=paste("TE: ", as.character(RM_cons[1,3]), "\n consensus size: ", as.character(RMcons_len), "bp; fragments: ", as.character(length(RM_cons$V1)), "; full length: ", as.character(length(full$V1))," (>=",as.character(as.numeric(FL_thresh)*100),"%)", sep = ""), cex.main = 0.9, xlab = "TE consensus (bp)", ylab = "divergence to consensus (%)")
    for(i in 1:length(RM_cons$V1)){
      segments(RM_cons$V5[i], RM_cons$V2[i], RM_cons$V6[i], RM_cons$V2[i], col=rgb(0,0,0,alpha=alpha))
    }
    for(i in 1:length(RM_cons$V1)){
      segments(full$V5[i], full$V2[i], full$V6[i], full$V2[i], col=rgb(1,0,0, alpha=full_alpha), lwd=1.5)
    }}
  else{
    par(mar=c(5,5,5,5))
    plot(range(0, RMcons_len), range(0, auto_y), type = "n", main=paste("TE: ", as.character(RM_cons[1,3]), "\n consensus size: ", as.character(RMcons_len), "bp; fragments: ", as.character(length(RM_cons$V1)), "; full length: ", as.character(length(full$V1))," (>=",as.character(as.numeric(FL_thresh)*100),"%)", sep = ""), cex.main = 0.9, xlab = "TE consensus (bp)", ylab = "divergence to consensus (%)")
    for(i in 1:length(RM_cons$V1)){
      segments(RM_cons$V5[i], RM_cons$V2[i], RM_cons$V6[i], RM_cons$V2[i], col=rgb(0,0,0,alpha=alpha))
    }
    for(i in 1:length(RM_cons$V1)){
      segments(full$V5[i], full$V2[i], full$V6[i], full$V2[i], col=rgb(1,0,0, alpha=full_alpha), lwd=1.5)
    }}
  
  #make the coverage matrix and graph
  cov=matrix(rep(0, length(RM_cons$V1)*as.numeric(RMcons_len)), byrow = T, ncol = as.numeric(RMcons_len))
  for(i in 1:length(RM_cons$V1)){
    times=RMcons_len-RM_cons$V6[i]
    cov[i,]<-c(rep(0,RM_cons$V5[i]-1),rep(1,abs(RM_cons$V5[i]-RM_cons$V6[i])+1), rep(0,times))
  }
  par(new=T)
  plot(colSums(cov), type='l', axes = F, ylab = NA, xlab = NA, col=covcol, ylim = c(0, max(colSums(cov))))
  axis(side = 4)
  mtext(side = 4, line = 3, 'consensus coverage (bp)')
}
