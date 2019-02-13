#' context_likelihood
#' 
#' Find the mutually similar cells between batches using the context likelihood
#' 
#' @param x a gene expression matrix (normalized and scaled)
#' @param batch a numeric vector indicating the batch labels of each cell
#' @param p0 cutoff for context likelihood neighbors (default: 0.2)
#'
#' @return a matrix of context likelihood neighbors (CLN)
#'
#' @author Wuming Gong
#'
context_likelihood <- function(x, batch, p0 = 0.2){

	M <- ncol(x)	# number of cells

	if (missing(batch))
		stop('batch is missing')

	batch <- as.numeric(factor(batch))
	bs <- unique(batch)	# unique batch indicators
	nbt <- length(bs)

	if (nbt == 1){
		D <- rdist(t(x))
	}else{
		D <- matrix(NA, M, M)	# distance matrix of cells at current time points
		param <- as.matrix(expand.grid(bs, bs))
		param <- param[param[, 1] != param[, 2], , drop = FALSE]
		for (i in 1:nrow(param)){
			m1 <- which(batch == param[i, 1])
			m2 <- which(batch == param[i, 2])
			D[m1, m2] <- rdist(t(x[, m1]), t(x[, m2]))
		}
	}

	# z-score of each cell's distance to all other non-batch cells
	Z <- (rowMeans(D, na.rm = TRUE) - D) / rowSds(D, na.rm = TRUE)	
	Z[Z < 0] <- 0
	Z[is.na(Z)] <- 0
	Z <- sqrt(Z^2 + t(Z^2)) # compute the pairwise z-score
	Z[Z < qnorm(1 - p0)] <- 0
	Z <- as(Z, 'dgCMatrix')
	Z

} # context_likelihood

