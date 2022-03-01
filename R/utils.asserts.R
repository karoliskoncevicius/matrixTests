assert_numeric_mat_or_vec <- function(x) {
  name <- as.character(substitute(x))
  if(is.null(x) || !is.numeric(x) | !(is.matrix(x) | is.vector(x)))
    stop('"', name, '"', ' must be a numeric matrix or vector')
}

assert_numeric_mat <- function(x) {
  name <- as.character(substitute(x))
  if(is.null(x) || !is.numeric(x) | !is.matrix(x))
    stop('"', name, '"', ' must be a numeric matrix')
}

assert_vec_length <- function(x, ...) {
  name   <- as.character(substitute(x))
  lens   <- unlist(list(...))
  lnames <- as.character(substitute(list(...)))[-1]
  lnames <- paste(lnames, collapse=' or ')
  if(!(length(x) %in% lens) | (NCOL(x) > 1 & NROW(x) > 1))
    stop('"', name, '"', ' must be a vector with length ', lnames)
}

assert_logical_vec_length <- function(x, ...) {
  name   <- as.character(substitute(x))
  lens   <- unlist(list(...))
  lnames <- as.character(substitute(list(...)))[-1]
  lnames <- paste(lnames, collapse=' or ')
  if(!(length(x) %in% lens) | !is.logical(x) | (NCOL(x) > 1 & NROW(x) > 1))
    stop('"', name, '"', ' must be a logical vector with length ', lnames)
}

assert_character_vec_length <- function(x, ...) {
  name   <- as.character(substitute(x))
  lens   <- unlist(list(...))
  lnames <- as.character(substitute(list(...)))[-1]
  lnames <- paste(lnames, collapse=' or ')
  if(!(length(x) %in% lens) | !is.character(x) | (NCOL(x) > 1 & NROW(x) > 1))
    stop('"', name, '"', ' must be a character vector with length ', lnames)
}

assert_numeric_vec_length <- function(x, ...) {
  name   <- as.character(substitute(x))
  lens   <- unlist(list(...))
  lnames <- as.character(substitute(list(...)))[-1]
  lnames <- paste(lnames, collapse=' or ')
  if(!(length(x) %in% lens) | !is.numeric(x) | (NCOL(x) > 1 & NROW(x) > 1))
    stop('"', name, '"', ' must be a numeric vector with length ', lnames)
}

assert_all_in_set <- function(x, vals) {
  name <- as.character(substitute(x))
  vnames <- paste(vals, collapse=", ")
  if(is.null(x) | !all(x %in% vals))
    stop('all "', name, '" values must be in: ', vnames)
}

assert_all_in_open_interval <- function(x, min, max) {
  name <- as.character(substitute(x))
  if(is.null(x) | any(anyNA(x) | x<=min | x>=max))
    stop('all "', name, '" values must be greater than ', min, ' and lower than ', max)
}

assert_all_in_closed_interval <- function(x, min, max) {
  name <- as.character(substitute(x))
  if(is.null(x) | any(anyNA(x) | x<min | x>max))
    stop('all "', name, '" values must be between: ', min, ' and ', max)
}

assert_equal_nrow <- function(x, y) {
  namex <- as.character(substitute(x))
  namey <- as.character(substitute(y))
  if(nrow(x) != nrow(y))
    stop('"', namex, '" and "', namey, '" must have the same number of rows')
}

assert_equal_ncol <- function(x, y) {
  namex <- as.character(substitute(x))
  namey <- as.character(substitute(y))
  if(ncol(x) != ncol(y))
    stop('"', namex, '" and "', namey, '" must have the same number of columns')
}

assert_ncol_equal_nrow <- function(x, y) {
  namex <- as.character(substitute(x))
  namey <- as.character(substitute(y))
  if(ncol(x) != nrow(y))
    stop('number of "', namex, '" columns must be equal to the number of "', namey, '" rows')
}

assert_max_number_of_levels <- function(x, mlevels) {
  name <- as.character(substitute(x))
  if(is.null(x) || length(stats::na.omit(unique(x))) > mlevels)
    stop('"', name, '"', ' must have no more than ', mlevels, ' unique elements')
}

assert_unique_colnames <- function(x) {
  name <- as.character(substitute(x))
  if(is.null(colnames(x)) || any(duplicated(colnames(x))))
    stop('column names of "', name, '"', ' must be unique')
}

assert_nested_model <- function(m0, m1) {
  if(!all(colnames(m0) %in% colnames(m1)) || !all.equal(m0, m1[,colnames(m0),drop=FALSE], check.attributes=FALSE)) {
    stop('null model must be nested (subset of the full model)')
  }
}
