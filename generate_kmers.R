generate_kmers <- function(dna, k) {
  len <- nchar(dna)
  starts <- 1:(len - k + 1)
  ends <- starts + k - 1
  kmers <- sapply(1:length(starts), function(i) substr(dna, starts[i], ends[i]))
  return(data.frame(
    seq = kmers,
    start = starts,
    end = ends
  ))
}
