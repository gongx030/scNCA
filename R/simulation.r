#' simulateTemporalRNAseq
#'
#' Simulating time series single cell RNA-seq data with batch effects
#' 
#' @param ncells a matrix of number of cells.  The number of columns represents the batches and the number of rows represents the lineages (H). 
#' @param time.points a matrix of distribution of batches along the time indices.  The number of columns represents the batches and the rowsrepresents time points(T).  
#' @param K the number of dimensions of the low dimensional biological subspace (default: 5).
#' @param batch.effect the strength of the added batch effect, which corresponds to the weight of the Gaussian random batch vector (default: 5).
#' @param N number of genes (default: 2000).
#
#' @return A SummarizedExperiment object of simulated temporal scRNA-seq data
#' 
#' @details The simulation of temporal scRNA-seq data follows three basic assumptions: 
#' (1) the developmental trajectories can be represented within a K-dimensional biological subspace.  
#' (2) We assumed that all H lineages arise from the same progenitors and differentiated toward different directions.  
#' Thus, to simulate the developmental trajectory of any lineage h∈{1,⋯,H} , we use the origin of the K-dimensional space as the starting point, randomly select another point in the K-dimensional space as the terminal point, and draw a lineage segment (a line segment) that connects the origin and the terminal point.  
#' This lineage segment therefore represents the developmental trajectory of lineage h on the K-dimensional biological subspace.  
#' (3) The developmental speed is constant.  Thus, the each lineage segment is evenly split into T parts, where each part is associate with the time index from 1 to T.  
#' To simulate the cells of lineage h from time t∈{1,⋯,T}, we randomly draw points from the lineage segments that correspond to lineage h from time t on the K-dimensional space.  Note that for balanced data, the points are drawn from the entire lineage segments, while for the imbalanced data, the points can only be drawn from the lineage segments covered by the corresponding time period.  Then to simulate high-dimensional gene expression, the K-dimensional points are projected to N-dimensional space by a randomly Gaussian matrix.  Similar to Haghverdi et al., the batch effects are incorporated by generating a Gaussian random vector for each batch and adding it to the gene expression matrix12.  
#'
#' @export
#'
#' @author Wuming Gong
#'
#' @examples
#' library(scNCA)
#' # Let us first simulate a balanced temporal scRNA-seq data with 2,000 genes, three lineages, 
#' # three batches and five time points. In the balanced data, each batch covers the cells from all lineages and time points. 
#' # We assume that there are 100 cells for each lineage/batch combination so there are in total 900 cells:
#' ncells <- cbind(
#' 	c(100, 100, 100),   # batch i
#'  c(100, 100, 100),   # batch ii
#'  c(100, 100, 100)    # batch iii
#' )
#' time.points <- cbind(
#' 	c(1, 1, 1, 1, 1),   # batch i
#' 	c(1, 1, 1, 1, 1),   # batch ii
#' 	c(1, 1, 1, 1, 1)    # batch iii
#' )
#' set.seed(1)
#' se <- simulateTemporalRNAseq(ncells = ncells, time.points = time.points)
#'
simulateTemporalRNAseq <- function(ncells, time.points, K = 5, batch.effect = 5, N = 2000){

	nb <- ncol(ncells)	# number of batches
	nc <- nrow(ncells)	# number of lineages
	M <- sum(ncells)
	nt <- nrow(time.points)	# number of time points
	time.points <- time.points %*% diag(1 / colSums(time.points))

	batch <- rep(1:nb, colSums(ncells))
	mem <- unlist(lapply(1:nb, function(b) rep(1:nc, ncells[, b])))

	C <- sparseMatrix(i = 1:M, j = mem, dims = c(M, nc))
	G <- sparseMatrix(i = 1:M, j = batch, dims = c(M, nb))

	# all clusters start from origin 
	# ends of each trajectory (cluster)
	ends <- matrix(rnorm(nc * K, sd = 5), ncol = K)	# nc by K 
	Mu <- matrix(NA, M, K)

	# determining the time index of each cell from each batch
	CT <- do.call('rbind', lapply(1:nb, function(b){
		m <- batch == b
		cbind(which(m), sample(1:nt, sum(m), time.points[, b], replace = TRUE))
	}))
	CT <- sparseMatrix(i = CT[, 1], j = CT[, 2], dims = c(M, nt))

	# randomly sample cells along the trajectory from the origin to the end
	for (c in 1:nc){	# for each cluster
		d0 <- sqrt(sum(ends[c, ]^2))
		b <- seq(0, 1, length.out = nt + 1)
		for (t in 1:nt){	# for each time point
			m <- mem == c & CT[, t]
			if (any(m)){
				ds <- runif(sum(m), min = d0 * (t - 1) / nt, max = d0 * t / nt)
				r <- ds / d0	# relative positions between start (origin) and end
				Mu[m, ] <- r %o% ends[c, ]
			}
		}
	}

	sigma <- matrix(rgamma(nc * K, 1, 1), ncol = K)	# nc by K
	V <- matrix(rnorm(M * K, mean = Mu, sd = sigma), M, K)

	proj <- matrix(rnorm(N * K), nrow = N, ncol = K)
	X <- proj %*% t(V)
	X <- X + rnorm(N * M)# Adding some random noise.

	# random gene level batch vectors
	W <- do.call('cbind', lapply(1:nb, function(b) batch.effect * rnorm(N)))
	X <- X + W %*% t(G)	# Adding a random batch vector

	se <- SummarizedExperiment(assays = list(counts = as.matrix(X)), colData = DataFrame(batch = batch, group = mem, time = max.col(CT)))
	colData(se)$time.table <- CT

	# Generate the the cell labels
	cn <- sprintf('%d_%d', batch, mem)
	colData(se)$cell_label <- cn

	se

} # simulateTemporalRNAseq
