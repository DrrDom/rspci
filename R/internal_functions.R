get_col <- function(x,y) {x[,if(is.name(substitute(y))) deparse(substitute(y)) else y, drop = FALSE][[1]]}



get_full_names <- function(df, addM = TRUE, addN = TRUE) {
  if (addM & addN) {
    return(dplyr::mutate(df[, c("FragID", "M", "N")], full_name = paste0(FragID, " (M = ", M, ", N = ", N, ")"))[, "full_name"])
  } else {
    if (addM) {
      return(dplyr::mutate(df[, c("FragID", "M")], full_name = paste0(FragID, " (M = ", M, ")"))["full_name"])
    } else {
      if (addN) {
        return(dplyr::mutate(df[, c("FragID", "N")], full_name = paste0(FragID, " (N = ", N, ")"))["full_name"])
      }
    }
  }
}



replace_na <- function(vec, value = 0) {
  vec[is.na(vec)] <- value
  return(vec)
}
