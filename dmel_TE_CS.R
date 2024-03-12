library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(tidyverse)
library(gridExtra)
# Load TE annotation in 18 Drosophila species
load("RData/edta18_te_anno.RData", verb = T)
# Load the function that makes a data frame from GRanges
source(".dffromgr.R")
# Load TE information table
load("/home/au732606/R/projects/STRIPE/RData/tes_and_hierarchy_3rd.RData",
     verb = T)
# Load Dmel genome
dmel_ont <- readDNAStringSet("/home/au732606/projects/DrosoTE/fifteen_genomes/Dmel/ont.wglkd.fasta")
# Function for getting conservation score across MSA consensus
getCS <- function(dnaset, mat){
  # MSA on supplied sequences
  ds_msa <- msa(dnaset, method = "ClustalW")
  # Make a consensus based on MSA
  consensus <- msaConsensusSequence(ds_msa,type = 'upperlower', ignoreGaps=T)
  # Calculate conservation score
  cs <- msaConservationScore(ds_msa, mat, gapVsGap=1)
  # Remove the gaps from consensus and CS coordinates
  cs_filt <- cs[gregexpr("[^-]", consensus)[[1]]]
  # Make a rolling window mean transformation of CS using 20 bp windows
  cs_sw20 <- rollapply(cs_filt, width = 20,
                   FUN = mean, by = 20, align = "left")
  cs_sw20.df <- cbind(seq(1, length(cs_filt), by = 20)[1:length(cs_sw20)],
               cs_sw20)
  colnames(cs_sw20.df) <- c("pos", "CS")
  return(cs_sw20.df)
}

# Matrix for conservation score
mat <- nucleotideSubstitutionMatrix(4, -1)
mat <- cbind(rbind(mat, "-"=-3), "-"=-3)

View(df_from_gr(te18$Dmel))
# Try CS calculation on one TE with manually approved insertions
# I won't curate all the insertions for general approach
mdg3 <- readDNAStringSet("fasta/dmel_mdg3.fa")
mdg3_msa <- msa(mdg3, method = "ClustalW")
cons_mdg3 <- msaConsensusSequence(mdg3_msa,type = 'upperlower', ignoreGaps=T)
mdg3_cs <- msaConservationScore(mdg3_msa, mat, gapVsGap=1)
plot(mdg3_cs[gregexpr("[^-]", cons_mdg3)[[1]]], type = 'h')
plot(mdg3_cs, type = 'h')
mdg3_cs_filt <- mdg3_cs[gregexpr("[^-]", cons_mdg3)[[1]]]
heh <- rollapply(mdg3_cs_filt, width = 20,
                 FUN = mean, by = 20, align = "left")
hoh <- cbind(seq(1, length(mdg3_cs_filt), by = 20)[1:length(heh)],
             heh)
colnames(hoh) <- c("pos", "CS")
p <- ggplot(as_tibble(hoh), aes(x = pos, y = CS/(4*length(mdg3)^2)))+
  geom_area()+
  theme_bw()+
  ylab("Conservation score")
ggsave(device = "pdf", "plots/mdg3_conserve_score_along_consensus.pdf", plot = p, width = 10, height = 4)

# Generate CS plots for all the TEs in our database
te.full.cs <- lapply(te3.hier$id, function(te){
  # Get a reference length
  te_length <- te3.hier$length[te3.hier$id == te]
  # Leave only insertions of specific TE that are close or are full-length
  te_gr <- te18$Dmel[te18$Dmel$Name == te & width(te18$Dmel) >= 0.9 * te_length &
                       width(te18$Dmel) <= 1.1 * te_length]
  # Let's build MSA only if we have at least 5 fl insertions
  if (length(te_gr) < 5) {
    return(NULL)
  } else {
    # Get sequences of insertions
    te_seq <- getSeq(dmel_ont, te_gr)
    print(te)
    # Calculate CS
    consscore <- getCS(te_seq, mat)
    # Prepare a plot, we divide CS by the max score to normalize: 1 means perfect identity
    p <- ggplot(as_tibble(consscore), aes(x = pos, y = CS/(4*length(te_seq)^2)))+
      geom_area()+
      theme_bw()+
      ylab("Conservation score")+
      ggtitle(paste0(te, " (", length(te_seq), " insertions)"))
    return(p)
  }
})

# Remove NULL elements from the list of plots
plot_list <- Filter(Negate(is.null), te.full.cs)
# Prepare multi page plot
multi_page_plot <- marrangeGrob(grobs = plot_list, nrow = 4, ncol = 1)

# Save to a PDF file
ggsave("plots/dmel_TE_CS.pdf", multi_page_plot, device = "pdf", width = 11, height = 8.5)  


