##' Create an empty sparse matrix over the field Z2
##'
##' The matrix is represented in a sparse format as a list of rows, where
##' each row is represented by as an integer vector storing the (sorted)
##' column numbers of the non-zero elements of the row.
##'
##' It is not necessary to specify the number of columns at the time the
##' matrix is created, it will be inferred from the maximum non-zero entry.
##'
##' @title Create Sparse Matrix over Z2
##' @param nrow Number of rows in the matrix
##' @param ncol Number of columns in the matrix
##' @return A `sparseZ2` matrix with all elements 0
##' @export
sparseZ2 <- function(nrow, ncol=NA) {
  structure(
    vector("list", nrow),
    Dim = as.integer(c(nrow, ncol)),
    class = "sparseZ2"
  )
}

##' Get the dimensions of a `sparseZ2` matrix
##'
##' @title Dimensions of a Sparse Matrix
##' @param x A `sparseZ2` matrix
##' @return A vector of length 2 giving the number of rows and columns.
##' @export
dim.sparseZ2 <- function(x) {
  attr(x, "Dim")
}

##' Print method for sparse matrices over Z2
##'
##' Print a `sparseZ2` matrix in standard matrix format
##'
##' @param x A `sparseZ2` matrix
##' @param ... Additional arguments passed to print
##' @export
print.sparseZ2 <- function(x, ...) {
  mat <- as.matrix(x)
  cat("Sparse matrix over Z2 (", nrow(x), "x", ncol(x), "):\n", sep="")
  print(mat, ...)
  invisible(x)
}

##' Set elements in a `sparseZ2` matrix.
##'
##' Set elements of a `sparseZ2` matrix, with any non-zero values treated as 1.
##' It is an error for the column index to be missing.
##'
##' @title Set Elements of a Sparse Matrix.
##' @param x A `sparseZ2` matrix
##' @param i Row index or vector of row indices
##' @param j Column index or vector of column indices
##' @param value Value(s) to set (0 or 1)
##' @return Updated sparse matrix
##' @export
`[<-.sparseZ2` <- function(x, i, j, value) {
  if(missing(j)) stop("Column index cannot be missing")
  if(missing(i)) i <- seq_len(length(x))
  if(length(i)==1) {
    x[[i]] <- as.integer(
      sort(unique(c(setdiff(x[[i]], j[value==0]), j[value!=0]))))
  } else {
    value <- matrix(value, length(i),length(j))
    for(k in seq_along(i))
      x[[i[k]]] <- as.integer(
          sort(unique(c(setdiff(x[[i[k]]], j[value[k,]==0]), j[value[k,]!=0]))))
  }
  x
}

##' Extract elements from a sparse matrix over the field Z2
##'
##' Extract elements from a `sparseZ2` matrix.  It is an error for the
##' column index to be missing.
##'
##' @title Extract Elements of a Sparse Matrix.
##' @param x A `sparseZ2` matrix
##' @param i Row index or vector of row indices
##' @param j Column index or vector of column indices
##' @return The value(s) at the specified position(s) (0 or 1)
##' @export
`[.sparseZ2` <- function(x, i, j) {
  if(missing(j)) stop("Column index cannot be missing")
  if(missing(i)) i <- seq_len(length(x))
  if(length(i)==1) {
    ifelse(j %in% x[[i]], 1, 0)
  } else {
    m <- matrix(0,length(i),length(j))
    for(k in seq_along(i)) m[k,] <- ifelse(j %in% x[[i[k]]], 1, 0)
    m
  }
}

##' Convert a matrix to a `sparseZ2` matrix.
##'
##' Converts a matrix to a `sparseZ2` format, with any non-zero values
##' treated as 1.
##'
##' @title Convert Matrix to Sparse Format.
##' @param mat A matrix of 0s and 1s
##' @return A sparse matrix representation
##' @export
as.sparseZ2 <- function(mat) {

  sparse <- sparseZ2(nrow(mat), ncol(mat))
  for (i in seq_len(nrow(mat))) {
    sparse[[i]] <- which(mat[i,] != 0)
  }
  sparse
}

##' Convert a sparse matrix to a regular matrix.
##'
##' @title Convert Sparse Matrix to Regular Matrix.
##' @param x A `sparseZ2` matrix
##' @param ... ignored.
##' @return A regular matrix of 0s and 1s
##' @export
as.matrix.sparseZ2 <- function(x,...) {
  dm <- dim(x)
  if(is.na(dm[2])) dm[2] <- do.call(max,x)
  mat <- matrix(0L, nrow = dm[1], ncol=dm[2])
  for (i in seq_along(x)) {
    mat[i, x[[i]]] <- 1L
  }
  mat
}

##' Multiply a `sparseZ2` matrix by a vector over the field Z2
##'
##' @title Sparse Matrix Vector Multiplication
##' @param A A `sparseZ2` matrix
##' @param b A vector over Z2
##' @return A vector over Z2
##' @export
multAb <- function(A,b) {
  vapply(A,function(r) as.integer(sum(b[r]) %% 2),0L)
}


##' Transpose a `sparseZ2` matrix
##'
##' @title Transpose a Sparse Matrix
##' @param A A `sparseZ2` matrix
##' @return The transposed `sparseZ2` matrix
##' @export
transpose <- function(A) {
  B <- sparseZ2(ncol(A), nrow(A))
  for(i in seq_len(nrow(B)))
    B[[i]] <- which(vapply(A, function(x) i %in% x, TRUE))
  B
}

##' Combine `sparseZ2` matrices by rows
##'
##' @title Combine `sparseZ2` Matrices by Rows
##' @param ... `sparseZ2` matrices to be combined
##' @param deparse.level Not used
##' @return A sparse matrix
##' @export
rbind.sparseZ2 <- function(..., deparse.level = 1) {
  args <- list(...)
  ncols <- vapply(args,ncol,0L)
  if(length(ncols <- ncols[!is.na(ncols)])==0) ncols <- NA
  A <- do.call(c,args)
  structure(A,
            Dim=as.integer(c(length(A),ncols))  ,
            class="sparseZ2")
}




##' Reduce two rows of a sparse matrix over the field Z2.
##'
##' Each row of the sparse matrix is represented by an integer vector
##' the elements of which are the (sorted) column numbers of the
##' non-zero elements of the row. This function takes sorted integer
##' vectors representing two rows of a matrix and computes the sorted
##' integer vector corresponding to sum of the rows in Z2.
##'
##' @title Row Reduction in Z2
##' @param r1 sorted integer vector
##' @param r2 sorted integer vector
##' @return a sorted integer vector representing the reduced row.
##' @useDynLib sparseZ2
##' @export
reduce_row <- function(r1,r2) {
  .Call("xor_merge",as.integer(r1),as.integer(r2))
}

#reduce_row <- function(r1, r2) {
#  as.integer(sort(union(setdiff(r1, r2), setdiff(r2, r1))))
#}




##' Calculates the reduced form of the sparse linear system A x = b by
##' Gaussian or Gauss-Jordan elimination over the field Z2.
##'
##' This implementation does no pivoting so the order of the rows and
##' columns is preserved.
##'
##' @title Sparse Gauss-Jordan Elimination over Z2
##' @param A A `sparseZ2` matrix
##' @param b A vector over the field Z2
##' @param jordan Whether Gauss-Jordan or Gaussian reduction is performed
##'
##' @return A list containing the reduced matrix and right hand side.
##'
##' @examples
##' # Create a sparse matrix
##' A <- sparseZ2(4, 4)
##' A[1,c(1,3,4)] <- 1
##' A[2,c(2,3)] <- 1
##' A[3,c(1,2,4)] <- 1
##' A[4,c(1,4)] <- 1
##'
##' b <- c(1, 0, 1, 1)
##'
##' # Perform Gaussian elimination
##' sys <- gaussJordanZ2(A, b, jordan = FALSE)
##' sys
##'
##' # Perform Gauss-Jordan elimination
##' sys <- gaussJordanZ2(A, b, jordan = TRUE)
##' sys
##'
##' @export
gaussJordanZ2 <- function(A, b, jordan = TRUE) {
  m <- length(A)

  ## Forward elimination
  for (i in seq_len(m)) {
    if(length(A[[i]])>0) {
      p <- A[[i]][1]
      if(jordan) {
        ## Reduce all rows above the pivot
        for(j in seq_len(i-1)) {
          if(p %in% A[[j]]) {
            A[[j]] <- reduce_row(A[[j]], A[[i]])
            b[j] <- (b[j] + b[i]) %% 2
          }
        }
      }

      ## Reduce all rows below the pivot
      for(j in seq_len(m-i)+i) {
        if(p %in% A[[j]]) {
          A[[j]] <- reduce_row(A[[j]],A[[i]])
          b[j] <- (b[j] + b[i]) %% 2
        }
      }

    }
  }
  list(A=A, b=as.integer(b))
}




##' Calculates the reduced form of the sparse linear system A x = b by
##' Gaussian elimination over the field Z2.
##'
##' This implmenetation does no pivoting so the order of the rows and
##' columns is preserved.
##'
##' @title Sparse Gaussian Elimination over Z2
##' @param A A sparse matrix over the field Z2 (see details)
##' @param b A vector over the field Z2
##'
##' @return A list containing the reduced matrix and right hand side.
##'
##' @examples
##' # Create a sparse matrix
##' A <- sparseZ2(4, 4)
##' A[1,c(1,3,4)] <- 1
##' A[2,c(2,3)] <- 1
##' A[3,c(1,2,4)] <- 1
##' A[4,c(1,4)] <- 1
##'
##' b <- c(1, 0, 1, 1)
##'
##' # Perform Gaussian elimination
##' sys <- gaussZ2(A, b)
##' sys
##'
##' @export
gaussZ2 <- function(A, b) {
  class(A) <- NULL
  m <- length(A)

  ## Leading column of a row
  lead_col <- function(x) {
    if(length(x)>0) x[1] else NA
  }

  ## Leading columns of the rows or NA if already used
  lead <- vapply(A, lead_col, 0L)

  ## Forward elimination
  while(any(!is.na(lead))) {
    i <- which.min(lead)
    p <- lead[i]
    lead[i] <- NA
    ks <- which(p==lead)
    A[ks] <- lapply(A[ks], reduce_row, A[[i]])
    b[ks] <- (b[ks] + b[i]) %% 2
    lead[ks] <- vapply(A[ks], lead_col, 0L)
  }
  class(A) <- "sparseZ2"
  list(A=A,b=as.integer(b))
}


##' Reduce a linear system in row echelon form to reduced row-echelon form.
##'
##' Given a Gaussian reduced system computed by [gaussJordanZ2()]
##' compute the full Gauss-Jordan reduction.
##'
##' @title Sparse Gauss-Jordan Elimination over Z2
##' @param A A sparse matrix over the field Z2
##' @param b A vector over the field Z2
##' @return A list containing the reduced matrix and right hand side
##' @examples
##' # Create a sparse matrix
##' A <- sparseZ2(4, 4)
##' A[1,c(1,3,4)] <- 1
##' A[2,c(2,3)] <- 1
##' A[3,c(1,2,4)] <- 1
##' A[4,c(1,4)] <- 1
##'
##' b <- c(1, 0, 1, 1)
##'
##' # Reduce to row echelon form by Gaussian elimination
##' sys <- gaussJordanZ2(A, b, jordan=FALSE)
##' sys
##'
##' # Reduce to reduced row echelon form
##' sys <- jordanReduceZ2(sys$A,sys$b)
##' sys
##'
##' # Reduce to reduced row echelon form by Gauss-Jordan elimination
##' sys <- gaussJordanZ2(A, b, jordan=TRUE)
##' sys
##'
##' @export
jordanReduceZ2 <- function(A, b) {
  m <- length(A)

  for(i in rev(seq_len(m))) {
    if(length(A[[i]]) > 0) {
      for(j in seq_len(i-1)) {
        if(A[[i]][1] %in% A[[j]]) {
          A[[j]] <- reduce_row(A[[j]], A[[i]])
          b[j] <- (b[j] + b[i]) %% 2
        }
      }
    }
  }

  list(A=A, b=b)
}

##' Calculates the reduced form of the sparse linear system A x = b by
##' Gaussian elimination over the field Z2.
##'
##' This implementation does no pivoting so the order of the rows and
##' columns is preserved.
##'
##' @title Sparse Gaussian Elimination over Z2
##' @param A A `sparseZ2` matrix
##' @param b A vector over the field Z2
##'
##' @return A list containing the reduced matrix and right hand side.
##'
##' @examples
##' # Create a sparse matrix
##' A <- sparseZ2(4, 4)
##' A[1,c(1,3,4)] <- 1
##' A[2,c(2,3)] <- 1
##' A[3,c(1,2,4)] <- 1
##' A[4,c(1,4)] <- 1
##'
##' b <- c(1, 0, 1, 1)
##'
##' # Perform Gaussian elimination
##' sys <- gaussZ2(A, b)
##' sys
##'
gaussZ2A <- function(A, b) {

  m <- length(A)
  b <- as.integer(b)

  ## Leading column of a row
  lead_col <- function(x) {
    if(length(x)>0) x[1] else NA_integer_
  }

  ## Leading columns of the rows or NA if already used
  lead <- vapply(A, lead_col, 0L)

  ## Forward elimination
  while(any(!is.na(lead))) {
    i <- which.min(lead)
    p <- lead[i]
    lead[i] <- NA
    for(j in which(p==lead)) {
      A[[j]] <- reduce_row(A[[j]], A[[i]])
      b[j] <- (b[j] + b[i]) %% 2
      lead[j] <- lead_col(A[[j]])
    }
  }

  list(A=A,b=as.integer(b))
}




##' Reduce a linear system using a previously reduced system
##'
##' Given a system A1 x1 = b1 that has been reduced with [gaussJordanZ2()] or
##' [gaussZ2()], reduce the system A2 x2 = b2 using the reduced form of A1.
##'
##' @title Sparse Gaussian Reduction over Z2
##' @param A1 The reduced sparse matrix from the first system
##' @param b1 The reduced right-hand side from the first system
##' @param A2 The sparse matrix to be reduced
##' @param b2 The right-hand side to be reduced
##'
##' @return A list containing the reduced A2 and b2
##' @examples
##' # Create a sparse matrix
##' A <- sparseZ2(4, 4)
##' A[1,c(1,3,4)] <- 1
##' A[2,c(2,3)] <- 1
##' A[3,c(1,2,4)] <- 1
##' A[4,c(1,4)] <- 1
##'
##' b <- c(1, 0, 1, 1)
##'
##' # Split the system
##' A1 <- sparseZ2(2, 4)
##' A1[1,c(1,3,4)] <- 1
##' A1[2,c(2,3)] <- 1
##' b1 <- b[1:2]
##'
##' A2 <- sparseZ2(2, 4)
##' A2[1,c(1,2,4)] <- 1
##' A2[2,c(1,4)] <- 1
##' b2 <- b[3:4]
##'
##' # Reduce A1
##' sys1 <- gaussZ2(A1, b1)
##' sys1
##'
##' # Reduce A2 using A1
##' sys2 <- gaussReduceZ2(sys1$A, sys1$b, A2, b2)
##' sys2
##'
##' # Further reduce A2
##' sys2 <- gaussZ2(sys2$A, sys2$b)
##' sys2
##'
##' @export
gaussReduceZ2 <- function(A1, b1, A2, b2) {
  m1 <- length(A1)
  m2 <- length(A2)
  n <- max(0, unlist(c(A1, A2)))

  ## Leading column of a row
  lead_col <- function(x) {
    if(length(x)>0) x[1] else NA_integer_
  }

  ## Leading columns of the rows or NA if already used
  lead <- vapply(A1, lead_col, 0L)

  ## Reduce A2 using pivots from A1
  for(i in order(lead)) {
    if(length(A1[[i]])>0) {
      for(j in seq_len(m2)) {
        if(A1[[i]][1] %in% A2[[j]]) {
          A2[[j]] <- reduce_row(A2[[j]], A1[[i]])
          b2[j] <- (b2[j] + b1[i]) %% 2
        }
      }
    }
  }

  list(A = A2, b = b2)
}



##' Find the solutions of a (possibly underdetermined) linear system over Z2.
##'
##' Given a reduced matrix A and vector b, `consistentSystem` checks if the
##' system is consistent and has at least one solution.  If the system is
##' consistent, `determinedVars` returns the variables that are determined by
##' the system, and `undeterminedVars` returns the variables that are
##' undetermined.  When `all=FALSE`, `undeterminedVars` returns only variables
##' that occur in the matrix.
##'
##' Given a Gauss-Jordan reduced matrix A and vector b that is consistent,
##' `solution` returns a one possible solution, with any undetermined variables
##' set to 0.
##'
##' @title Solutions of Linear Systems
##' @param A A reduced sparse matrix over Z2
##' @param b A reduced vector over Z2
##' @param all Include all variables or only those that appear in the matrix.
##' @return TRUE if the system is consistent, FALSE otherwise
##' @examples
##' # Create a sparse matrix
##' A <- sparseZ2(4, 4)
##' A[1,c(1,3,4)] <- 1
##' A[2,c(2,3)] <- 1
##' A[3,c(1,2,4)] <- 1
##' A[4,c(1,4)] <- 1
##'
##' b <- c(1, 0, 1, 1)
##'
##' ## Determine solutions
##' sys <- gaussJordanZ2(A, b, jordan = FALSE)
##' consistentSystem(sys$A,sys$B)
##' undeterminedVars(sys$A)
##'
##' sys <- jordanReduceZ2(sys$A, sys$b)
##' solution(sys$A, sys$b)
##'
##' @rdname solution
##' @export
solution <- function(A, b) {
  b <- as.integer(b)
  ks <-vapply(unclass(A)[lengths(A)>0 & b!=0], function(r) r[1], 0L)
  n <- ncol(A)
  if(is.na(n)) n <- max(unlist(A))
  x <- integer(n)
  x[ks] <- b[ks]
  x
}

##' @rdname solution
##' @export
consistentSystem <- function(A, b) {
  all(b[lengths(A)==0]==0)
}

##' @rdname solution
##' @export
determinedVars <- function(A) {
  sort(vapply(unclass(A)[lengths(A)>0], function(r) r[1], 0L))
}

##' @rdname solution
##' @export
undeterminedVars <- function(A,all=FALSE) {
  A1 <- unclass(A)[lengths(A)>0]
  sols <- vapply(A1, function(r) r[1], 0L)
  if(all) {
    n <- ncol(A)
    if(is.na(n)) n <- max(unlist(A))
    vars <- setdiff(seq_len(n), sols)
  } else {
    vars <- unique(unlist(lapply(A1, function(r) r[-1])))
  }
  sort(setdiff(vars, sols))
}


##' Compute the null space of a sparse matrix over Z2
##'
##' Given a reduced matrix A and vector b that are consistent, `nullSpace`
##' returns a sparse matrix with rows that are the non-trivial elements of a
##' basis for the null space of A.
##'
##' @title Null Space of a Sparse Matrix
##' @param A A sparse matrix over Z2
##' @return A sparse matrix with rows that form a basis for the null space of A
##' @examples
##' # Create a sparse matrix
##' A <- sparseZ2(4, 4)
##' A[1,c(1,3,4)] <- 1
##' A[2,c(2,3)] <- 1
##' A[3,c(1,2,4)] <- 1
##' A[4,c(1,4)] <- 1
##'
##' b <- c(1, 0, 1, 1)
##'
##' ## Find solutions
##' sys <- gaussJordanZ2(A, b, jordan = FALSE)
##' consistentSystem(sys$A,sys$B)
##' solution(sys$A,sys$b)
##' nullSpace(sys$A)
##' @export
nullSpace <- function(A) {
  vs <- undeterminedVars(A,all=TRUE)
  B <- lapply(vs, function(v) {
    ks <- vapply(A, function(r) if(v %in% r) r[1] else NA, 0L)
    sort(c(v,ks[!is.na(ks)]))
  })
  structure(B,Dim=c(length(B),NA),class="sparseZ2")
}
