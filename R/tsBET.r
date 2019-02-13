#' tsBET
#'
#' Testing for batch effets in temporal scRNA-seq data
#' batch effect test for temporal scRNA-seq data
#' 
#' @param x a M by M sparse matrix representing the adjacency graph of the cells
#' @param batch a numeric vector indicating the batch labels of each cell
#' @param time.table a binary matrix indicating the distribution of cells at each time point
#' @param resample number of resampling for building the null distribution of local batch score (default: 100)
#' @param neighborhood_size number of neighboring cells when computing the local batch score (default: 100)
#' @param mc.cores number of CPU cores (default: 8)
#'
#' @return a list of tsBET results
#'
#' @export
#' 
#' @author Wuming Gong
#' 
tsBET <- function(x, batch, time.table, resample = 100, neighborhood_size = 100, mc.cores = 8){

	if (missing(batch))
		stop('batch must be specified')

	batch <- as.numeric(factor(batch))
	nb <- max(batch)	# number of batches

	M <- nrow(x)	# number of cells

	if (missing(time.table))
		time.table <- Matrix(TRUE, M, 1)

	nt <- ncol(time.table)	# number of time point
	time <- max.col(time.table)

	B <- sparseMatrix(i = 1:M, j = batch, dims = c(M, nb))	# cell ~ batch

	g <- graph.adjacency(x, mode = 'undirected', weighted = TRUE)
	g_mst <- mst(g)	# the minimum spanning tree

	BT <- as(t(B), 'dgCMatrix') %*% time.table	# batch ~ time

	# the cells not farther than a given limit from give cells
	# aka. the neighborhood of each cell
	neighborhood <- lapply(ego(g_mst, neighborhood_size), as.numeric)

	res <- do.call('rbind', mclapply(1:length(neighborhood), function(i){
		m <- neighborhood[[i]]	# neighboood cells of cell m
		Mi <- length(m)
		batchh <- batch[m]	# the batch index
		timeh <- time[m]	# the time index
		ah <- summary(x[m, m])	# the sub-tree

		# randomly sampling batch index based on the time index within the neighborhood
		U <- matrix(NA, Mi, resample)
		for (t in 1:nt){
			mt <- timeh == t
			if (any(mt)){
				U[mt, ] <- sample(1:nb, sum(mt) * resample, replace = TRUE, prob = BT[, t])
			}
		}

		# number of adjacency cell pairs that are from the same batch in the real data
		s <- sum(batchh[ah[, 1]] == batchh[ah[, 2]]) / nrow(ah)	

		# number of adjacent cell pairs that are from the same batch from the resampling data
		s0 <- colSums(U[ah[, 1], ] == U[ah[, 2], ]) / nrow(ah)

		data.frame(
			cell = i,
			neighborhood_size = Mi,
			ds = s - s0
		)
	}, mc.cores = mc.cores))

	pval <- 1 - (sum(res[, 'ds'] > 0) + 1) / (nrow(res) + 1)
	list(p.value = pval, neighborhood_size = neighborhood_size, resample = resample, M = M)
	
} # tsBET



