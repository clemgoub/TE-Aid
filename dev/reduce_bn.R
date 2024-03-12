beppo <- read.table("/home/au732606/projects/DrosoTE/fifteen_genomes/Dsuz/Circe/circe_cand.fa.genome.blastn.out",
                    colClasses = c("character", "character", "numeric", rep("integer", 7), "numeric", "numeric")) 
colnames(beppo) <- c("query", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "score")
beppo <- beppo %>% arrange(
  subject, sstart, -send, qstart
)
Rcpp::sourceCpp(file.path("~/Software/TE-Aid", "reduce.cpp"))

reduce_blast <- function(blast_df){
  bstart <- blast_df$sstart
  bdf <- blast_df
  bdf$sstart[bstart > bdf$send] <- bdf$send[bstart > bdf$send]
  bdf$send[bstart > bdf$send] <- bstart[bstart > bdf$send]
  bdf <- bdf[order(bdf$subject,
                        bdf$sstart,
                        -bdf$send,
                        bdf$qstart), ]
  nored <- lapply(split(bdf, bdf$subject), function(df){
    if (nrow(df) == 1){
      return(df)
    } else {
      index <- sapply(2:nrow(df), function(idx){
        df$send[idx] >= df$send[idx - 1]
      })
      return(df[index,])
    }
  })
  return(do.call(rbind, nored))
}
blast_red <- reduce_blastn_results(beppo)
blast_slowred <- reduce_blast(beppo)
all(beppo$send <= beppo$sstart)
