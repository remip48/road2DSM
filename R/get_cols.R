get_cols <- function(mat, stem, K) mat[, paste0(stem, "[", seq_len(K), "]"), drop = FALSE]
