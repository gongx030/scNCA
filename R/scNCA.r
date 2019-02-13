#' scNCA
#'
#' @import Matrix
#' @import SummarizedExperiment
#' @importFrom parallel mclapply 
#' @importFrom irlba irlba prcomp_irlba partial_eigen ssvd
#' @importFrom matrixStats rowSds rowVars rowMedians
#' @importFrom fields rdist
#' @importFrom FNN knnx.index knn.index
#' @importFrom gplots colorpanel
#' @importFrom cluster pam
#' @importFrom igraph graph.adjacency mst ego
#' @importFrom S4Vectors DataFrame
#' @import tensorflow
#' @import futile.logger
#' @name scNCA
NULL
# > NULL


#' scNCA
#' 
#' Integrating multiple sources of temporal scRNA-seq data by neighborhood component analysis
#' 
#' @param x a gene expression matrix (normalized and scaled)
#' @param K number of low dimensions
#' @param batch a numeric vector indicating the batch labels of each cell
#' @param time.table a binary matrix indicating the distribution of cells at each time point
#' @param max.iter maximum iterations (default: 5000)
#' @param p0 cutoff for context likelihood neighbors (default: 0.2)
#' @param learning_rate learing rate for AdaGrad (default: 0.01)
#'
#' @return a list of scNCA results. 
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
#' X <- assays(se)$counts	# the simulated scRNA-seq expression matrix
#' X <- as.matrix(X %*% Diagonal(x = 1 / sqrt(colSums(X^2))))	# consine scaling
#'
#' set.seed(1)
#' mf <- scNCA(X, batch = colData(se)$batch, time.table = colData(se)$time.table)
#' 
#' # set up the color and shape of the plot
#' n2bg <- c(
#'  '1_1' = 'red', '1_2' = 'lightblue', '1_3' = 'green',  # batch i
#'  '2_1' = 'red', '2_2' = 'lightblue', '2_3' = 'green',  # batch ii
#'  '3_1' = 'black', '3_2' = 'black', '3_3' = 'black'     # batch iii
#' )
#' n2pch <- c(
#'  '1_1' = 21, '1_2' = 21, '1_3' = 21, # batch i
#'  '2_1' = 24, '2_2' = 24, '2_3' = 24, # batch ii
#'  '3_1' = 3, '3_2' = 3, '3_3' = 3     # batch iii
#' )
#' n2col <- c(
#'  '1_1' = 'black', '1_2' = 'black', '1_3' = 'black',  # batch i
#'  '2_1' = 'black', '2_2' = 'black', '2_3' = 'black',  # batch ii
#'  '3_1' = 'red', '3_2' = 'lightblue', '3_3' = 'green' # batch iii
#' )
#' bg <- n2bg[colData(se)$cell_label]
#' pch <- n2pch[colData(se)$cell_label]
#' col2 <- n2col[colData(se)$cell_label]
#' 
#' library(DDRTree)
#' dt <- DDRTree(mf$Y, verbose = TRUE, sigma = 1e-3, maxIter = 5)
#' y <- t(dt$Z)
#' plot(y[, 1], y[, 2], bg = bg, pch = pch, col = col2, cex = 1.5, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
#'
scNCA <- function(x, K = 5, batch = NULL, time.table = NULL, max.iter = 5000, p0 = 0.20, learning_rate = 0.01){

	N <- nrow(x)
	M <- ncol(x)	# number of cells

	if (is.null(batch))
		batch <- rep(1, M)
	else
		batch <- as.numeric(factor(batch))


	if (is.null(time.table))
		time.table <- as(matrix(TRUE, nrow = M, ncol = 1), 'ngCMatrix')

	time.table <- time.table[, colSums(time.table) > 0, drop = FALSE]	# remove time points w/o any cells
	time <- max.col(time.table)
	nt <- ncol(time.table)	# number of time points
	nb <- max(batch)	# number of batches

	# the transformation matrix for each batch
	A <- lapply(1:nb, function(h) matrix(runif(K * N, min = -1, max = 1), K, N))
	A <- lapply(1:nb, function(h) tf$Variable(tf$cast(A[[h]], tf$float32)))

	loss <- 0
	for (t in 1:nt){	# for each time point
		m <- time == t	# cells at current time point
		Mt <- sum(m)	# number of total cells at current time point
		batcht <- batch[m]	# batch labels for cells at current time point

		# reorder the cells based on the batch IDs
		# so that it is easier for tensorflow to manipulate data
		co <- order(batcht)
		xt <- x[, m]
		xt <- xt[, co]
		batcht <- batcht[co]

		# compute the context likelihood neighbors
		W <- as.matrix(context_likelihood(xt, batch = batcht, p0 = p0))	
		W <- tf$cast(W, tf$float32)

		Yt <- lapply(1:nb, function(h){
			mh <- batcht == h
			if (any(mh))
				tf$matmul(A[[h]], tf$cast(xt[, which(mh)], tf$float32))
		})
		Yt <- Yt[!sapply(Yt, is.null)]	# remove non-existing batch at current time
		Y <- tf$concat(Yt, axis = 1L)

		# euclidean distance of each cell pairs on the transformed space
		YY <- tf$tile(tf$reduce_sum(Y^2, axis = 0L, keep_dims = TRUE), list(Mt, 1L))
		D <- YY + tf$transpose(YY) - 2 * tf$matmul(Y, Y, transpose_a = TRUE)

		# softmax over euclidean distance
		loss <- loss + tf$reduce_sum(-W * tf$nn$softmax(-D))
	}

	optimizer <- tf$train$AdagradOptimizer(learning_rate)
	train <- optimizer$minimize(loss)

	sess <- tf$Session()
	sess$run(tf$global_variables_initializer())

	for (iter in 1:max.iter){
		sess$run(train)
		if (iter == 1 || iter %% 500 == 0){
			cat(sprintf('[%s] iter=%5.d | loss=%10.3f\n', Sys.time(), iter, sess$run(loss)))
		}
	}

	A <- sess$run(A)
	Y <- matrix(0, K, M)
	for (h in 1:nb){
		m <- batch == h
		Y[, m] <- A[[h]] %*% x[, m]
	}

	list(Y = Y)

} # scNCA


