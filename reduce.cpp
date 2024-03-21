#include <Rcpp.h>
using namespace Rcpp;

// Helper function to swap start and end if start > end
void normalize_start_end(int &start, int &end) {
  if (start > end) {
    std::swap(start, end);
  }
}

// [[Rcpp::export]]
DataFrame reduce_blastn_results(DataFrame blastn) {
  CharacterVector subjects = blastn["subject"];
  IntegerVector sstarts = blastn["sstart"];
  IntegerVector sends = blastn["send"];
  IntegerVector qstarts = blastn["qstart"];
  // Additional columns
  CharacterVector queries = blastn["query"];
  NumericVector pidents = blastn["pident"];
  IntegerVector lengths = blastn["length"];
  IntegerVector mismatch = blastn["mismatch"];
  IntegerVector gapopen = blastn["gapopen"];
  IntegerVector qends = blastn["qend"];
  NumericVector evalues = blastn["evalue"];
  NumericVector scores = blastn["score"];
  
  int n = subjects.size();
  
  // Normalize sstart and send
  for (int i = 0; i < n; ++i) {
    normalize_start_end(sstarts[i], sends[i]);
  }
  
  // Sort based on subject, sstart, and sends
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int i, int j) {
    if (subjects[i] < subjects[j]) return true;
    if (subjects[i] > subjects[j]) return false;
    if (sstarts[i] < sstarts[j]) return true;
    if (sstarts[i] > sstarts[j]) return false;
    return sends[i] > sends[j];  // send in descending order
  });
  
  // Merge and remove redundant ranges
  std::vector<int> filtered_idx;
  for (int i = 0; i < n; ++i) {
    int current_idx = idx[i];
    bool include = true;
    
    // Compare with previous entries
    for (int j : filtered_idx) {
      if (subjects[j] == subjects[current_idx] && 
          sstarts[j] <= sstarts[current_idx] && 
          sends[j] >= sends[current_idx]) {
        // Current range is within a previous range
        if (sstarts[j] < sstarts[current_idx] || sends[j] > sends[current_idx]) {
          include = false;
          break;
        } else if (qstarts[j] <= qstarts[current_idx]) {
          // Current range is identical to a previous range but has a higher qstart
          include = false;
          break;
        }
      }
    }
    
    if (include) {
      filtered_idx.push_back(current_idx);
    }
  }
  
  // Create vectors for the resulting DataFrame
  CharacterVector result_queries(filtered_idx.size());
  CharacterVector result_subjects(filtered_idx.size());
  NumericVector result_pidents(filtered_idx.size());
  IntegerVector result_lengths(filtered_idx.size());
  IntegerVector result_mismatch(filtered_idx.size());
  IntegerVector result_gapopen(filtered_idx.size());
  IntegerVector result_qstarts(filtered_idx.size());
  IntegerVector result_qends(filtered_idx.size());
  IntegerVector result_sstarts(filtered_idx.size());
  IntegerVector result_sends(filtered_idx.size());
  NumericVector result_evalues(filtered_idx.size());
  NumericVector result_scores(filtered_idx.size());
  
  for (size_t i = 0; i < filtered_idx.size(); ++i) {
    result_queries[i] = queries[filtered_idx[i]];
    result_subjects[i] = subjects[filtered_idx[i]];
    result_pidents[i] = pidents[filtered_idx[i]];
    result_lengths[i] = lengths[filtered_idx[i]];
    result_mismatch[i] = mismatch[filtered_idx[i]];
    result_gapopen[i] = gapopen[filtered_idx[i]];
    result_qstarts[i] = qstarts[filtered_idx[i]];
    result_qends[i] = qends[filtered_idx[i]];
    result_sstarts[i] = sstarts[filtered_idx[i]];
    result_sends[i] = sends[filtered_idx[i]];
    result_evalues[i] = evalues[filtered_idx[i]];
    result_scores[i] = scores[filtered_idx[i]];
  }
  
  // Create the resulting DataFrame
  return DataFrame::create(
    Named("query") = result_queries,
    Named("subject") = result_subjects,
    Named("pident") = result_pidents,
    Named("length") = result_lengths,
    Named("mismatch") = result_mismatch,
    Named("gapopen") = result_gapopen,
    Named("qstart") = result_qstarts,
    Named("qend") = result_qends,
    Named("sstart") = result_sstarts,
    Named("send") = result_sends,
    Named("evalue") = result_evalues,
    Named("score") = result_scores
  );
}
