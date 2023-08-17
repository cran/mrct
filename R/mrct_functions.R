#' Pairwise inner product for \eqn{L^2} functions
#'
#' Calculate all pairwise inner products between elements from \eqn{L^2} supplied to this function. The integral is approximated by the Trapezoidal rule for uniform grids:
#' \deqn{ \int_a^b f(x) dx \approx \Delta x \left( \sum_{i=1}^{N-1} f(x_i) + \frac{f(x_N) - f(x_0)}{2} \right) }
#' whereas \eqn{ \{x_i \}} is an uniform grid on \eqn{[a,b]} such that \eqn{a = x_0 < x_1 < \ldots < x_N = b} and \eqn{\Delta x} the step size, i.e. \eqn{\Delta x := x_2 - x_1}.
#' Therefore, it is assumed that the functions are evaluated at the same, equidistant grid.
#' @param grid A numeric vector of the uniform grid on which the functions are evaluated.
#' @param data A numeric matrix. Each function has to be a vector stored in a column of `data` and evaluated at the points of `grid`.
#' Thus, the number of rows and columns of `data` correspond to `length(grid)` and the number of functions, respectively.
#' @return Numeric symmetric matrix containing the approximated pairwise inner products between the functions supplied by `data`. The entry \eqn{(i,j)} of the result is the inner product
#' between the \eqn{i}-th and \eqn{j}-th column of `data`.
#' @examples
#' # Create orthogonal fourier basis via `fdapace` package
#' library(fdapace)
#' basis <- fdapace::CreateBasis(K = 10,
#'                               type = "fourier")
#' iP <- innerProduct(grid = seq(0, 1, length.out = 50), # default grid in CreateBasis()
#'                    data = basis)
#' round(iP,3)
#' # Since the basis is orthogonal, the resulting matrix will be the identity matrix.
#' @export
innerProduct <- function(grid,
                         data
){

  step <- diff(grid[1:2])
  nfunc <- dim(data)[2]
  grid.length <- length(grid)

  M <- matrix(0,nfunc,nfunc)

  for (i in 1:nfunc){
    for (j in 1:nfunc){
      func_prod <- data[,i] * data[,j]
      M[i,j] <- step * ( sum(func_prod[2:(grid.length-1)]) + (func_prod[1]+func_prod[grid.length])/2 )
    }
  }
  return(M)
}

mrct.rescale <- function(data,
                         data.cov,
                         alpha,
                         H,
                         scaling.iterations = 5,
                         scaling.tolerance = 10^(-4)
                         ){

  k0 <- iter <- 1
  k1 <- Inf
  error <- (k0-k1)^2

  eigen.cov <- eigen(data.cov)

  data.sqrtcov <- eigen.cov$vectors %*% diag(abs(eigen.cov$values)^(1/2)) %*% t(eigen.cov$vectors)

  std.data <- t(solve(data.cov + alpha/k0*diag(rep(1,ncol(data)))) %*%  data.sqrtcov %*% t(data))
  std.data.cov <- std.data %*% t(std.data)
  mhd <- diag(std.data.cov)

  while(iter <= scaling.iterations & error > scaling.tolerance){
    std.data <- t(solve(data.cov + alpha/k0*diag(rep(1,ncol(data)))) %*%  data.sqrtcov %*% t(data))
    std.data.cov <- std.data %*% t(std.data)
    mhd <- diag(std.data.cov)

    k1 <- stats::median(mhd)/stats::median(mrct.distmhd(data.cov=data.cov,
                                     alpha=alpha,
                                     k=k0))
    error <- (k1-k0)^2
    iter <- iter + 1
    k0 <- k1
  }

  return(list("aMHD" = mhd,
              "std.data" = std.data,
              "hsubset" = order(mhd)[1:H],
              "scalingparameter" = k1)
         )
}

mrct.distmhd <- function(data.cov,
                         alpha,
                         k = 1,
                         iter = 2000,
                         seed =123){
  set.seed(seed)

  p <- dim(data.cov)[1]
  Q <- data.cov %*% solve(data.cov + alpha/k*diag(rep(1,p)))
  eigv.Q <- eigen(Q,symmetric = T)$values

  realizations <- c()
  for(i in 1:iter){
    beta <- stats::rchisq(p,1)
    realizations[i] <- sum(eigv.Q^2 * beta)
  }
  return(realizations)
}

#' Minimum regularized covariance trace estimator
#'
#' Functional outlier detection based on the `minimum` `regularized` `covariance` `trace` estimator \insertCite{oguamalam2023minimum}{mrct} as a robust covariance estimator.
#' This estimator uses a generalization of the Mahalanobis distance for the functional setting \insertCite{berrendero2020}{mrct} and a corresponding theoretical cutoff value.
#' @references \insertRef{berrendero2020}{mrct}.
#' @references \insertRef{oguamalam2023minimum}{mrct}.
#' @importFrom Rdpack reprompt
#' @param data Numeric matrix of a functional data set for which the esimator has to be calculated. Each row contains an observation. They are assumed to be observed on the same regular grid.
#' @param h Numeric value between \eqn{0.5} and \eqn{1}. Ratio of the data which the estimator is based on. Default is set to \eqn{0.75}, i.e. \eqn{75\%} of the data will be used for the estimator.
#' @param alpha Numeric (default is \eqn{0.01}). Tikhonov regularization parameter \eqn{\alpha}.
#' @param initializations Integer (default is \eqn{5}). Number of random initial subsets.
#' @param subset.iteration Integer (default is \eqn{10}). Maximum number of how often each subset is re-estimated and adjusted.
#' @param seed Integer (default is \eqn{123}). Random seed for reproducibility.
#' @param scaling.iterations Integer (default is \eqn{5}). The maximum number of times \eqn{k_1} is re-scaled if the error between subsequent
#' scalingparameters does not fall below `scaling.tolerance`.
#' @param scaling.tolerance Numeric (default is \eqn{10^{-4}}). The error tolerance for re-scaling. If the error falls below this value, the re-scaling procedure stops.
#' @param criterion Character. Criterion based on which the optimal subset is chosen among the final subsets. Possible options are: "`cluster`" and the default "`sum`".
#' @param sum.percentage Numeric value between \eqn{0.5} and \eqn{1}. Corresponding to the "`sum`" criterion. Determines the fraction of observations up to which the sum over the sorted functional Mahalanobis distances is calculated (in ascending order). Default is set to \eqn{0.75}, i.e. the sum of the smallest \eqn{75\%} of Mahalanobis distances is calculated. If outliers are present, this value should not be to high, in order not to include any outlying curves.
#' @return A list:
#' \item{theoretical}{Integer vector of the indices corresponding to the outliers based on the MRCT estimator.}
#' \item{theoretical.w}{Same as `theoretical` with an additional re-weighting step.}
#' \item{aMHD}{Numeric vector containing the functional Mahalanobis distances of all observations based on the MRCT estimator.}
#' \item{aMHD.w}{Same as `aMHD` with an additional re-weighting step.}
#' \item{quant}{Numeric. Theoretical cutoff value for outlier detection.}
#' \item{quant.w}{Same as `quant` with an additional re-weighting step.}
#' \item{k}{Numeric. Scalingparameter \eqn{k_1} of Algorithm 1 described in \insertCite{oguamalam2023minimum}{mrct}.}
#' \item{k.w}{Same as `k` with an additional re-weighting step.}
#' \item{optimal.subset}{Integer vector of the optimal h-subset.}
#' \item{subsets}{Numeric matrix containing all final subsets. Each row of `subsets` is one final subset.}
#' \item{objval}{Numeric vector with the objective values of the final subsets based on `criterion`.}
#' @examples
#' # Fix seed for reproducibility
#' set.seed(123)
#'
#' # Sample outlying indices
#' cont.ind <- sample(1:50, size=10)
#'
#' # Generate 50 curves on the interval [0,1] at 50 timepoints with 20% outliers
#' y <- mrct.rgauss(x.grid=seq(0,1,length.out=50), N=50, model=1,
#'                  outliers=cont.ind, method="linear")
#'
#' # Visualize curves (regular curves grey, outliers black)
#' colormap <- rep("grey",50); colormap[cont.ind] <- "black"
#' matplot(x=seq(0,1,length.out=50), y=t(y), type="l", lty="solid",
#'         col=colormap, xlab="t",ylab="")
#'
#' # Run MRCT
#' mrct.y <- mrct(data=y, h=0.75, alpha=0.1,
#'                initializations=10, criterion="sum")
#'
#' # Visualize alpha-Mahalanobis distance with cutoff (horizontal black line)
#' # Colors correspond to simulated outliers, shapes to estimated (MRCT) ones
#' # (circle regular and triangle irregular curves)
#' shapemap <- rep(1,50); shapemap[mrct.y$theoretical.w] <- 2
#' plot(x=1:50, y=mrct.y$aMHD.w, col=colormap, pch=shapemap,
#'      xlab="Index", ylab=expression(alpha*"-MHD"))
#' abline(h = mrct.y$quant.w)
#'
#' # If you dont have any information on possible outliers,
#' # alternatively you could use the S3 method plot.mrct()
#' mrct.plot(mrct.y)
#' @export
mrct <- function(data,
                 h = 0.75,
                 alpha = 0.01,
                 initializations = 5,
                 subset.iteration = 10,
                 seed = 123,
                 scaling.iterations = 10,
                 scaling.tolerance = 10^(-4),
                 criterion = "sum",
                 sum.percentage = 0.75
){
  if(h < 0 | h > 1){
    stop("The value of `h` cannot be negative or above 1!")
  }
  if(h < 0.5){
    stop("`h` = ",h," was selected. For values below 0.5 the results might break down if too many outliers are present. Please consider a value between 0.5 and 1.")
  }
  if(sum.percentage < 0 | sum.percentage > 1){
    stop("The value of `sum.percentage` cannot be negative or above 1!")
  }

  N <- nrow(data)
  p <- ncol(data)

  # Allocation of data frames
  objval <- c()
  hsubsets <- matrix(NA,
                     ncol=floor(h*N),
                     nrow=initializations)
  scalingparameters <- c()

  ######################################################################################
  ### Part 1: Start with many initial subsets and iterate until convergence for each ###
  ######################################################################################
  for(i in 1:initializations){

    set.seed(seed+i)
    subset.initial <- sample(1:N,floor(h*N))

    data.centered <- data - matrix(rep(robustbase::colMedians(data),N), ncol = p, byrow = T)
    data.cov <- t(data.centered[subset.initial,]) %*% data.centered[subset.initial,] / length(subset.initial)

    k <- 1

    subset.old <- subset.initial
    subset.new <- Inf

    while(k <= subset.iteration
          & setequal(subset.old,subset.new) == F
          ){

      if(k >= 2){subset.old <- subset.new}

      # Calculate new h-subset
      tmp <- mrct.rescale(data = data.centered,
                          data.cov = data.cov,
                          H = floor(N*h),
                          alpha = alpha,
                          scaling.iterations = scaling.iterations
                          )

      subset.new <- tmp$hsubset

      data.centered <- data - matrix(rep(colMeans(data[subset.new,]),N),byrow=T,ncol=p)
      data.cov <- stats::cov(data.centered[subset.new,])

      k <- k + 1
    }

    mhd <- tmp$aMHD / tmp$scalingparameter

    if(criterion == "sum"){
      cutoff <- stats::quantile(mhd,sum.percentage)
      objval[i] <- sum(mhd[mhd <= cutoff])
    }else if(criterion == "cluster"){
      res <- stats::kmeans(mhd,centers = 2)
      objval[i] <- res$tot.withinss/res$betweenss
    }

    hsubsets[i,] <- subset.new
    scalingparameters[i] <- tmp$scalingparameter
  }

  ###################################
  ### Part 2: Select final subset ###
  ###################################
  # Take subset with smallest obj value
  subset.optimal <- hsubsets[order(objval)[1],]
  scalingparameter.optimal <- scalingparameters[order(objval)[1]]

  # "Extract" optimal covariance and scaled a-MHD
  data.centered <- data - matrix(rep(colMeans(data[subset.optimal,]),N),byrow=T,ncol=p)
  data.cov <- stats::cov(data.centered[subset.optimal,])
  eigen.cov <- eigen(data.cov)
  data.sqrtcov <- eigen.cov$vectors %*% diag(abs(eigen.cov$values)^(1/2)) %*% t(eigen.cov$vectors)
  std.data <- t(solve(data.cov + alpha/scalingparameter.optimal*diag(rep(1,p))) %*% data.sqrtcov %*% t(data.centered))
  mhd <- diag(std.data %*% t(std.data)) / scalingparameter.optimal

  dist <- mrct.distmhd(data.cov=data.cov,
                       alpha = alpha,
                       k = scalingparameter.optimal)

  quant <- stats::quantile(dist,0.99) # quantile

  ##############################
  ### Part 3: Weighting step ###
  ##############################
  non.outliers <- as.numeric(which(mhd <= quant))

  data.centered.weighted <- data - matrix(rep(colMeans(data[non.outliers,]),N),byrow=T,ncol=p)
  data.cov.weighted <- stats::cov(data.centered.weighted[non.outliers,])

  temp <- mrct.rescale(data = data.centered.weighted,
                       data.cov = data.cov.weighted,
                       alpha = alpha,
                       H = floor(h*N),
                       scaling.tolerance = scaling.tolerance,
                       scaling.iterations = scaling.iterations
                       )

  mhd.weighted <- temp$aMHD/temp$scalingparameter
  dist.weighted <- mrct.distmhd(data.cov = data.cov.weighted,
                                alpha = alpha,
                                k = temp$scalingparameter)
  quant.weighted <- stats::quantile(dist.weighted,0.99)
  # End of weighting

  ###########################################################
  ### Part 4: Determine outliers (with/without weighting) ###
  ###########################################################
  theoretical.outliers.weighted <- as.numeric(which(mhd.weighted >  quant.weighted))
  theoretical.outliers <- as.numeric(which(mhd >  quant))

  output <- list("theoretical" = theoretical.outliers,
                 "theoretical.w" = theoretical.outliers.weighted,
                 "aMHD" = mhd,
                 "aMHD.w" = mhd.weighted,
                 "quant" = quant,
                 "quant.w" = quant.weighted,
                 "k" = scalingparameter.optimal,
                 "k.w" = temp$scalingparameter,
                 "optimal.subset" = subset.optimal,
                 "subsets" = hsubsets,
                 "objval" = objval
  )

  return(output)
}

#' Sparse minimum regularized covariance trace estimator
#'
#' Robust outlier detection for sparse functional data as a generalization of the `minimum` `regularized` `covariance` `trace` (MRCT) estimator \insertCite{oguamalam2023minimum}{mrct}. At first the observations are smoothed
#' by a B-spline basis and afterwards the MRCT algorithm is performed with the matrix of basis coefficients.
#' @references \insertRef{oguamalam2023minimum}{mrct}.
#' @importFrom Rdpack reprompt
#' @inheritParams mrct
#' @param data Numeric matrix of a functional data set for which the esimator has to be calculated. Each row contains an observation. They are assumed to be observed on the same (probably sparse) regular grid. The number of grid points must be at least `nbasis`.
#' @param nbasis Integer. Number of B-spline basis functions for smoothing. The basis will be of order \eqn{4} and therefore, cannot contain less than \eqn{4} functions. The default value will be set to `dim(data)[2]`. i.e. the number of time points with a maximum of \eqn{15}.
#' @param new.p Integer. Length of the grid of the smoothed curves. The resulting grid will be an equidistant partition of `[rangeval[1],rangeval[length(rangeval)]]`. Default value is `dim(data)[2]`
#' @return A list with two entries
#' \item{mrct.output}{List. The same output as the function [mrct::mrct()]. For more details, see there.}
#' \item{data.smooth}{Numeric matrix. Collection of the smoothed curves of `data` with `dim(data)[1]` rows and `new.p` columns. Each row corresponds to one observation.}
#' @examples
#' # Fix seed for reproducibility
#' set.seed(123)
#'
#' # Sample outlying indices
#' cont.ind <- sample(1:50,size=10)
#'
#' # Generate 50 sparse curves on the interval [0,1] at 10 timepoints with 20% outliers
#' y <- mrct.rgauss(x.grid=seq(0,1,length.out=10), N=50, model=1,
#'                  outliers=cont.ind, method="linear")
#'
#' # Visualize curves (regular curves grey, outliers black)
#' colormap <- rep("grey",50); colormap[cont.ind] <- "black"
#' matplot(x = seq(0,1,length.out=10), y = t(y), type="l", lty="solid",
#'         col=colormap, xlab="t",ylab="")
#'
#' # Run sparse MRCT
#' sparse.mrct.y <- mrct.sparse(data = y, nbasis = 10, h = 0.75, new.p = 50,
#'                              alpha = 0.1, initializations = 10, criterion = "sum" )
#'
#' # Visualize smoothed functions
#' matplot(x=seq(0,1,length.out=50), y=t(sparse.mrct.y$data.smooth),
#'         type="l", lty="solid", col=colormap, xlab="t", ylab="")
#'
#' # Visualize alpha-Mahalanobis distance with cutoff (horizontal black line)
#' # Colors correspond to simulated outliers, shapes to estimated (sparse MRCT) ones
#' # (circle regular and triangle irregular curves)
#' shapemap <- rep(1,50); shapemap[sparse.mrct.y$mrct.output$theoretical.w] <- 2
#' plot(x = 1:50, y = sparse.mrct.y$mrct.output$aMHD.w, col=colormap, pch = shapemap,
#'      xlab = "Index", ylab = expression(alpha*"-MHD"))
#' abline(h = sparse.mrct.y$mrct.output$quant.w)
#'
#' # If you dont have any information on possible outliers,
#' # alternatively you could use the S3 method plot.mrctsparse()
#' mrct.sparse.plot(mrct.sparse.object = sparse.mrct.y)
#' @export
mrct.sparse <- function(data,
                        nbasis = dim(data)[2],
                        new.p = dim(data)[2],
                        h = 0.75,
                        alpha = 0.01,
                        initializations = 5,
                        seed = 123,
                        scaling.iterations = 10,
                        scaling.tolerance = 10^(-4),
                        criterion = "sum",
                        sum.percentage = 0.75
){

  N <- nrow(data)
  p <- ncol(data)

  rangeval <- 1:dim(data)[2]
  length.range <- length(rangeval)

  ##################################################################################
  #### Part 1: Smooth the data and extract the corresponding basis coefficients ####
  ##################################################################################

  # Smoothing is performed on B-spline basis function
  # Basis as well as matrix of coefficients are orthogonalized "by hand" after construction

  nbasis <- ifelse(nbasis > 15,15,nbasis)

  basis.coeffs <- matrix(0,
                         ncol = nbasis,
                         nrow = N)

  basis <- fda::create.bspline.basis(rangeval = c(rangeval[1],
                                                  rangeval[length.range]),
                                     nbasis = nbasis)

  for(i in 1:N){
    curve <- data[i,]
    curve.smooth <- fda::smooth.basis(argvals = c(rangeval[1]:rangeval[length.range])[!is.na(curve)],
                                      y = curve[!is.na(curve)],
                                      fdParobj = basis,
                                      covariates = NULL,
                                      method = "chol")
    basis.coeffs[i,] <- curve.smooth$fd$coefs
  }
  basis.coeffs[basis.coeffs < 0] <- 0
  data.smooth <- fda::fd(t(basis.coeffs),
                    basisobj = basis)
  data.smooth <- t(fda::eval.fd(evalarg = seq(rangeval[1],
                                              rangeval[length.range],
                                              length.out=new.p),
                                fdobj = data.smooth))

  ##################################################
  #### Part 2: Orthogonalize basis coefficients ####
  ##################################################

  # Evaluate basis on finer grid
  x.grid <- seq(from = rangeval[1],
                to = rangeval[length.range],
                length.out = 1000)
  basis.finer <- fda::eval.basis(evalarg = x.grid,
                                 basisobj = basis)
  # Matrix of inner products M
  M <- innerProduct(x.grid,
                          basis.finer)
  # Square root of M
  eigen.M <- eigen(M)
  M.sqrt <- eigen.M$vectors %*% diag(abs(eigen.M$values)^(1/2)) %*% t(eigen.M$vectors)

  # Orthogonal basis and corresponding coefficients
  basis.orthogonal <- solve(M.sqrt) %*% t(fda::eval.basis(evalarg = seq(rangeval[1],
                                                                   rangeval[2],
                                                                   length.out=p),
                                                     basisobj = basis))
  basis.coeffs.orthogonal <- basis.coeffs %*% M.sqrt

  ##########################################################
  #### Part 3: Run MRCT based on matrix of coefficients ####
  ##########################################################

  mrct.output <- mrct(data = basis.coeffs.orthogonal,
                      h = h,
                      alpha = alpha,
                      initializations = initializations,
                      seed = seed,
                      scaling.iterations = scaling.iterations,
                      scaling.tolerance = scaling.tolerance,
                      criterion = criterion,
                      sum.percentage = sum.percentage
  )
  output <- list("mrct.output" = mrct.output,
                 "data.smooth" = data.smooth)
  return(output)
}

mrct.kernel <- function(d,
                        sigma = 1,
                        l = 1,
                        method = "gaussian"){

  if(method == "gaussian"){
    return(sigma^2*exp(-d^2/(2*l^2)))
  }else if(method == "quadratic"){
    return(sigma*exp(-d^2/l))
  }else if(method == "linear"){
    return(sigma*exp(-abs(d)/l))
  }

}

#' Random sample from Gaussian process
#'
#' Generate random samples of Gaussian process on a uniform grid for different settings of the simulation study in \insertCite{oguamalam2023minimum;nobrackets}{mrct}.
#' @references \insertRef{oguamalam2023minimum}{mrct}.
#' @importFrom Rdpack reprompt
#' @param x.grid Numeric vector containing a uniform grid on which the process is defined.
#' @param N Integer number of observations to generate.
#' @param seed Integer (default is \eqn{123}).. Random seed for reprocudibility.
#' @param model Integer. Either \eqn{1, 2} or \eqn{3}. Corresponds to one of the three simulation settings.
#' @param outliers Integer vector containing the indices of outliers. If empty, then only regular curves will be generated.
#' @param sigma,l Numeric values with default equal to \eqn{1}. Parameters for the covariance kernel.
#' @param method Different types of covariance kernels. Possible options are "`quadratic`"
#' \deqn{\gamma(s,t) = \sigma \text{exp}\left(\frac{-(s-t)^2}{l}\right),}
#' "`linear`"
#' \deqn{\gamma(s,t) = \sigma \text{exp}\left(\frac{-|s-t|}{l}\right)}
#' or "`gaussian`" (default)
#' \deqn{\gamma(s,t) = \sigma^2 \text{exp}\left(\frac{-(s-t)^2}{2 l^2}\right)}.
#' @return Numeric matrix with \eqn{N} rows and `length(x.grid)` columns containing the randomly generated curves following a Gaussian process.
#' Each observations is a row of the result.
#' @examples
#' # Fix seed for reproducibility
#' set.seed(123)
#'
#' # Sample outlying indices
#' cont.ind <- sample(1:50,size=10)
#'
#' # Generate 50 curves on the interval [0,1] at 50 timepoints with 20% outliers
#' y <- mrct.rgauss(x.grid=seq(0,1,length.out=50), N=50 ,model=1,
#'                  outliers=cont.ind)
#'
#' # Visualize curves (regular curves grey, outliers black)
#' colormap <- rep("grey",50); colormap[cont.ind] <- "black"
#' matplot(x=seq(0,1,length.out=50), y=t(y), type="l", lty="solid",
#'         col=colormap, xlab="t",ylab="")
#' @export
mrct.rgauss <- function(x.grid,
                        N,
                        seed=123,
                        model,
                        outliers,
                        sigma = 1,
                        l = 1,
                        method = "linear"
                        ){
  set.seed(seed)
  grid.length <- length(x.grid)
  outliers.count <- length(outliers)

  distmat <- stats::dist(x.grid,diag=T,upper=T)
  distmat <- as.matrix(distmat)

  # Generate regular data
  K <- mrct.kernel(d = distmat,
                   sigma = sigma,
                   l = l,
                   method = method )

  L <- chol(K)

  e <- matrix(stats::rnorm(N*grid.length),grid.length,N)

  if(model == 1){
    main <- 30*x.grid*(1-x.grid)^1.5
  }else{
    main <- 4*x.grid
  }

  Y <- matrix(main,
              nrow=N,
              ncol=grid.length,
              byrow=T) + t(e)%*%L

  # Generate outlying curves
  if(outliers.count > 0){
    e <- matrix(stats::rnorm(outliers.count*grid.length),grid.length,outliers.count)

    if(model == 1){
      main <- 30*x.grid^1.5*(1-x.grid)

      Y[outliers,] <- matrix(main,nrow=outliers.count,ncol=grid.length,byrow=T) + t(e)%*%L

    }else if(model ==2){
      u <- as.numeric(stats::runif(outliers.count) < .5)
      mu <- stats::runif(outliers.count,.25,.75)
      t1 <- matrix(4*x.grid,
                   ncol=grid.length,
                   nrow=outliers.count,
                   byrow = T)
      t2 <- matrix(x.grid,
                   ncol=grid.length,
                   nrow=outliers.count,
                   byrow = T) - mu
      main <- t1 +  (-1)^u*1.8 + (0.02*pi)^(-1/2)*exp(-t2^2/0.02)

      Y[outliers,] <- main + t(e)%*%L

    }else{
      mu <- stats::runif(outliers.count,.25,.75)
      t1 <- matrix(4*x.grid,
                   ncol=grid.length,
                   nrow=outliers.count,
                   byrow = T)
      t2 <- matrix(x.grid,
                   ncol=grid.length,
                   nrow=outliers.count,
                   byrow = T) + mu
      main <- t1 + 2*sin(4*t2*pi)

      Y[outliers,] <- main + t(e)%*%L
    }

  }

  return(Y)
}

#' Integrated square error
#'
#' Calculates the approximation of the integrated square error between the estimated covariance based
#' on non-outlying curves of a data set determined by the MRCT estimator and the true kernel for one of the three outlier settings in the simulation study of \insertCite{oguamalam2023minimum;nobrackets}{mrct}.
#'
#' @references \insertRef{oguamalam2023minimum}{mrct}.
#' @importFrom Rdpack reprompt
#' @inheritParams mrct
#' @param outliers.est Integer vector containing the indices of outliers.
#' @param model Integer. \eqn{1} correspond to the first outlier setting, whereas \eqn{2} and \eqn{3} are related to the remaining two, which both have the same kernel.
#' @return Numeric value containing the approximated integrated square error between estimated and theoretical covariance.
#' @examples
#' # Fix seed for reproducibility
#' set.seed(124)
#'
#' # Sample outlying indices
#' cont.ind <- sample(1:100,size=10)
#'
#' # Generate 100 curves on the interval [0,1] at 150 timepoints with 20% outliers.
#' y <- mrct.rgauss(x.grid=seq(0,1,length.out=150), N=100, model=1,
#'                  outliers=cont.ind, method="linear")
#' # Run MRCT
#' mrct.y <- mrct(data=y, h=0.75, alpha=0.1,
#'                initializations=10, criterion="sum")
#' # Two additional curves are regarded as outlying according to the algorithm
#' mrct.y$theoretical.w %in% cont.ind
#' # Compare the ISE between true kernel and 1) true non-outliers,
#' # 2) estimated non-outliers and 3) the complete data
#' ise1 <- mrct.ise(data=y, outliers.est=cont.ind, model=1)
#' ise2 <- mrct.ise(data=y, outliers.est=mrct.y$theoretical.w, model=1)
#' ise3 <- mrct.ise(data=y, outliers.est=c(), model=1)
#' ise1; ise2; ise3
#'
#' @export
mrct.ise <- function(data,
                     outliers.est,
                     model){

  if(length(outliers.est) == 0){
    cov.est <- stats::cov(data)
  }else{
    cov.est <- stats::cov(data[-outliers.est,])
  }

  p <- dim(data)[2]
  x.grid <- seq(0,1,length.out = p)
  diff.grid <- as.matrix(stats::dist(x.grid,diag = T,upper = T))

  if(model == 1){
    cov.true <- 0.3*exp(-abs(diff.grid)/0.3)
  }else{
    cov.true <- exp(-abs(diff.grid))
  }

  return(sum((cov.est-cov.true)^2/p^2))
}


#' Plot function for result from [mrct::mrct()]
#'
#' A function for descriptive plots for an object resulting from a call to [mrct::mrct()].
#'
#' @param mrct.object A result from a call to [mrct::mrct()].
#' @importFrom fdapace CreateBasis
#' @return Descriptive plots
#' \item{aMHD.plot}{Alpha-Mahalanobis distances, corresponding cutoff values and coloring according to estimated outliers (grey regular, black irregular).}
#' \item{aMHD.plot.w}{Same as `aMHD.plot`, with additional re-weighting step.}
#' @examples
#' # Similar to example in mrct() helpfile
#' # Fix seed for reproducibility
#' set.seed(123)
#'
#' # Sample outlying indices
#' cont.ind <- sample(1:50, size=10)
#'
#' # Generate 50 curves on the interval [0,1] at 50 timepoints with 20% outliers
#' y <- mrct.rgauss(x.grid=seq(0,1,length.out=50), N=50, model=1,
#'                  outliers=cont.ind, method="linear")
#'
#' # Visualize curves (regular curves grey, outliers black)
#' colormap <- rep("grey",50); colormap[cont.ind] <- "black"
#' matplot(x=seq(0,1,length.out=50), y=t(y), type="l", lty="solid",
#'         col=colormap, xlab="t",ylab="")
#'
#' # Run MRCT
#' mrct.y <- mrct(data=y, h=0.75, alpha=0.1,
#'                initializations=10, criterion="sum")
#'
#' # Visualize alpha-Mahalanobis distance
#' # Colorinfromation according to estimated outliers (grey regular, black irregular)
#' mrct.plot(mrct.y)
#' @export
mrct.plot <- function(mrct.object){

  aMHD <- aMHD.w <- col.w <- value <- NULL

  # Plot of the alpha-Mahalanobis distances
  colormap <- colormap.w <- rep("grey",length(mrct.object$aMHD))
  colormap[mrct.object$theoretical] <- "black"
  colormap.w[mrct.object$theoretical.w] <- "black"

  data <- data.frame("aMHD" = mrct.object$aMHD,
                     "col" = colormap,
                     "aMHD.w" = mrct.object$aMHD.w,
                     "col.w" = colormap.w)

  plot1 <- ggplot2::ggplot(data = data, ggplot2::aes(x = 1:dim(data)[1],
                                    y = aMHD)) +
    ggplot2::geom_point(ggplot2::aes(col=col)) +
    ggplot2::scale_color_manual(values = c("black","grey")) +
    ggplot2::labs(x = "Index",
                  y = expression(alpha*"-MHD")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_hline(yintercept = mrct.object$quant)

  plot2 <- ggplot2::ggplot(data = data, ggplot2::aes(x = 1:dim(data)[1],
                                            y = aMHD.w)) +
    ggplot2::geom_point(ggplot2::aes(col=col.w)) +
    ggplot2::scale_color_manual(values = c("black","grey")) +
    ggplot2::labs(x = "Index",
                  y = expression("weighted "*alpha*"-MHD")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_hline(yintercept = mrct.object$quant.w)

  return(list("aMHD.plot"=plot1,
              "aMHD.plot.w"=plot2))
}

#' Plot function for result from [mrct::mrct.sparse()]
#'
#' A function for descriptive plots for an object resulting from a call to [mrct::mrct.sparse()].
#'
#' @param mrct.sparse.object A result from a call to [mrct::mrct.sparse()].
#' @param x Gridpoints on which the smoothed data is to be plotted on. Default is `seq(0,1,length.out=new.p)` whereas `new.p` is a parameter set in the call to [mrct::mrct.sparse()].
#' @return Descriptive plots.
#' \item{aMHD.plot}{Alpha-Mahalanobis distances, corresponding cutoff values and coloring according to estimated outliers (grey regular, black irregular).}
#' \item{aMHD.plot.w}{Same as `aMHD.plot`, with additional re-weighting step.}
#' \item{data.plot}{Plot of the smoothed curves, colors corresponding to estimated outliers (gery regular, black irregular). Per default, the x-axis is plotted over
#' an equidistant grid of the interval \eqn{[0,1]}. }
#' @examples
#' # Similar to example in mrct.sparse() helpfile
#' # Fix seed for reproducibility
#' set.seed(123)
#'
#' # Sample outlying indices
#' cont.ind <- sample(1:50,size=10)
#'
#' # Generate 50 sparse curves on the interval [0,1] at 10 timepoints with 20% outliers
#' y <- mrct.rgauss(x.grid=seq(0,1,length.out=10), N=50, model=1,
#'                  outliers=cont.ind, method="linear")
#'
#' # Visualize curves (regular curves grey, outliers black)
#' colormap <- rep("grey",50); colormap[cont.ind] <- "black"
#' matplot(x = seq(0,1,length.out=10), y = t(y), type="l", lty="solid",
#'         col=colormap, xlab="t",ylab="")
#'
#' # Run sparse MRCT
#' sparse.mrct.y <- mrct.sparse(data = y, nbasis = 10, h = 0.75, new.p = 50,
#'                              alpha = 0.1, initializations = 10, criterion = "sum" )
#'
#' # Visualize alpha-Mahalanobis distances and smoothed curves
#' # Colorinformation according to estimated outliers (grey regular, black irregular)
#' mrct.sparse.plot(mrct.sparse.object = sparse.mrct.y)
#'
#' @export
mrct.sparse.plot <- function(x = seq(0,1,length.out=dim(mrct.sparse.object[[2]])[2]), mrct.sparse.object){

  aMHD <- aMHD.w <- col.w <- Var1 <- Var2 <- value <- NULL

  mrct.object <- mrct.sparse.object[[1]]
  # Plot of the alpha-Mahalanobis distances
  colormap <- colormap.w <- rep("grey",length(mrct.object$aMHD))
  colormap[mrct.object$theoretical] <- "black"
  colormap.w[mrct.object$theoretical.w] <- "black"

  data <- data.frame("aMHD" = mrct.object$aMHD,
                     "col" = colormap,
                     "aMHD.w" = mrct.object$aMHD.w,
                     "col.w" = colormap.w)

  plot1 <- ggplot2::ggplot(data = data, ggplot2::aes(x = 1:dim(data)[1],
                                                     y = aMHD)) +
    ggplot2::geom_point(ggplot2::aes(col=col)) +
    ggplot2::scale_color_manual(values = c("black","grey")) +
    ggplot2::labs(x = "Index",
                  y = expression(alpha*"-MHD")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_hline(yintercept = mrct.object$quant)

  plot2 <- ggplot2::ggplot(data = data, ggplot2::aes(x = 1:dim(data)[1],
                                                     y = aMHD.w)) +
    ggplot2::geom_point(ggplot2::aes(col=col.w)) +
    ggplot2::scale_color_manual(values = c("black","grey")) +
    ggplot2::labs(x = "Index",
                  y = expression("weighted "*alpha*"-MHD")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_hline(yintercept = mrct.object$quant.w)

  # Plot smoothed curves
  # Colors correspond to estimated outlier information (grey regular, black irregular)
  curves <- mrct.sparse.object[[2]]
  data.c <- reshape2::melt(curves)
  x.cut <- floor(length(x)/6)
  x.label <- round(c(x[1],x[1+x.cut],x[1+2*x.cut],x[1+3*x.cut],x[+4*x.cut],x[length(x)]),2)
  plot3 <- ggplot2::ggplot(data = data.c, ggplot2::aes(x = Var2,
                                                       y = value,
                                                       fill = Var1)) +
    ggplot2::geom_line(col = "grey") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "t",
                  y = "") +
    ggplot2::scale_x_continuous(labels = x.label)

  return(list("aMHD.plot" = plot1,
              "aMHD.plot.w" = plot2,
              "data.plot" = plot3))
}
