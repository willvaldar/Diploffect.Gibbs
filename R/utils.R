
make.dip2hap.matrix <- function(H, prob.format, hap.names=as.character(1:H), sep="."){
  # check args
  stopifnot(is.wholenumber(H) & H > 1)
  stopifnot(length(hap.names)==H)
  stopifnot(prob.format %in% c("upper.tri", "lower.tri"))
  # first make multiplier matrix in lower.tri format
  mat <- matrix(0, nrow=H*(H+1)/2, ncol=H)
  dip.names <- character(nrow(mat))
  # diagonals
  diag(mat) <- 2
  dip.names[1:H] <- paste(hap.names, sep=sep, hap.names)
  # off-diagonals: construct table
  dips <- matrix(0, nrow=nrow(mat)-H, ncol=2)
  i <- 1
  for (h in 1:(H-1)){
    for (j in (h+1):H){
      dips[i, ] <- c(h, j)
      i <- i+1
    }
  }
  for (i in 1:nrow(dips)){
    mat[H+i, dips[i, ] ] <- 1
    dip.names[H+i] <- paste(collapse=sep, hap.names[dips[i, ]])
  }
  rownames(mat) <- dip.names
  colnames(mat) <- hap.names
  if ("upper.tri"==prob.format){
    ii <- reorder.diplotypes(H, to="upper.tri")
    mat <- mat[ii, ]
  } else if ("lower.tri"!=prob.format){
    stop("Unrecognized format '", prob.format, "'\n")
  }
  mat
}

#' Reorder a vector of diplotype probabilities from upper triangle format to lower triangle format.
#'
#' @description
#' For J=8 haplotypes, eg, A, B, C, ..., H, there are J*J=64 possible phased haplotype pairs (ie, diplotypes), which if written out in an J*J matrix would have rows (AA, AB, ..., AH), (BA, BB, ..., BH), etc. In typical applications considered here, however, diplotypes are unphased such that AB and BA are indistinguishable. ...
#' such that only the diagonal (AA, BB, ..., HH) and either the upper triangle or the lower triangl
#' 
reorder.diplotypes <- function(J, to.format, return.table=FALSE){
  # make a dummy matrix of all diplotype combinations
  str.mat <- matrix("", nrow=J, ncol=J)
  for (j in 1:J){
    for (k in 1:J){
      str.mat[j, k] <- paste(collapse="_", sort(c(j, k)))
    }
  }
  # collect the two common orderings
  upper <- c(diag(str.mat), str.mat[upper.tri(str.mat, diag=FALSE)])
  lower <- c(diag(str.mat), str.mat[lower.tri(str.mat, diag=FALSE)])
  tab <- list(upper=upper, lower=lower)
  if (return.table){
    return (tab)
  }
  # convert between the two orderings
  if ("lower.tri"==to.format){
    ii <- match(tab$lower, tab$upper)
  } else if ("upper.tri"==to.format){
    ii <- match(tab$upper, tab$lower)
  } else {
    stop(" 'to.format' argument must be 'lower.tri' or 'upper.tri' not '", to.format, "'")
  }
  ii
}
