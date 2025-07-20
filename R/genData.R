#' @title Generation of simulation data
#' @description \code{genData} This function generates simulation data (v0.1.0)
#'
#' @importFrom expm sqrtm
#' @importFrom dplyr arrange
#' @importFrom magrittr %>% divide_by extract multiply_by_matrix set_names
#'     set_colnames use_series
#' @importFrom purrr array_tree exec map map_dbl
#' @importFrom stats rnorm runif
#'
#' @param n0 the common sample size for each group
#' @param p the number of explanatory variables including intercept
#' @param m the number of groups (one of 10, 20 and 50)
#' @param case true cluster pattern (1 or 2)
#' @param rho a parameter of autocorrelation matrix for a matrix of expnalatory
#'   variables
#' @param SNR SNR value

#' @return a list with the following elements:
#' \item{y}{a vector of a response variable}
#'
#' \item{X}{a matrix of explanatory variables without intercept}
#'
#' \item{group}{a vector of group indexes}
#'
#' \item{adj}{a matrix with two columns expressing adjacency among groups}
#'
#' \item{Mu}{a vector of true mean}
#'
#' \item{XI}{a matrix of coefficients for cluster}
#'
#' \item{BETA}{a matrix of coefficients for group}
#'
#' \item{cluster}{a list expressing group indexes of each cluster}
#'
#' \item{cluster.labels}{a vector of cluster indexes}
#'
#' \item{SNR}{SNR}
#'
#' \item{SDC}{a dataframe with pairs of groups that they reside different clusters
#'     and standardized difference of coefficients
#' }
#'
#' @export
#' @examples
#' #genData(n0, p)

genData <- function(n0, p, m=10, case=1, rho=0.5, SNR=NULL){

  n <- n0*m
  sig <- 1; sig2 <- sig^2

  #---   group  ----------------------------------------------------------------

  if(m == 10)
  {
    m. <- ifelse(case==1, 3, 6)
  } else if(m == 20)
  {
    m. <- ifelse(case==1, 6, 12)
  } else if(m == 50)
  {
    m. <- ifelse(case==1, 10, 30)
  }

  TC <- true.cluster(m, m.)

  adj <- map(1:m, ~cbind(.x, setdiff(1:m, .x))) %>% exec(rbind, !!!.) %>%
    set_colnames(c("group", "adjacent"))

  group <- rep(1:m, each=n0)

  cluster <- map(1:m., ~cbind(c=.x, g=TC[[.x]])) %>%
    exec(rbind, !!!.) %>% as.data.frame %>% arrange(g) %>%
    use_series(c)

  E <- cbind(cluster[adj[,1]], cluster[adj[,2]]) %>% apply(1, sort) %>%
    t %>% unique %>% magrittr::extract(.[,1] - .[,2] != 0,) %>%
    array_tree(1)

  #---   mean structure   ------------------------------------------------------

  Psisq <- corM(p-1, rho)$sqrtM
  X <- runif(n*(p-1), -1, 1) %>% matrix(n, p-1) %>%
    magrittr::multiply_by_matrix(Psisq) %>% cbind(1, .)

  XI <- mapply(rep, 1:m., p) %>% t

  rownames(XI) <- paste0("c", 1:m.)

  BETA <- XI[cluster, ]

  X_ <- as.data.frame(X) %>% split(group)

  SNR0 <- map_dbl(E, ~{
    crossprod(Psisq%*%(XI[.x[1], -1] - XI[.x[2], -1])) %>%
      drop %>% divide_by(3*(p-1)*sig2)
  }) %>% mean

  if(is.null(SNR))
  {
    SNR <- SNR0
  } else
  {
    v <- sqrt(SNR / SNR0)
    XI <- XI*v; BETA <- BETA*v

    SNR <- map_dbl(E, ~{
      crossprod(Psisq%*%(XI[.x[1], -1] - XI[.x[2], -1])) %>%
        drop %>% divide_by(3*(p-1)*sig2)
    }) %>% mean
  }

  Mu <- map(1:m, ~{
    data.matrix(X_[[.x]])%*%BETA[.x,] %>% drop
  }) %>% unlist

  y <- Mu + rnorm(n, sd=sig)

  idx <- sample(1:n, n)

  return(list(
    y = y[idx], X = X[idx,-1], group = group[idx], adj = adj,
    Mu = Mu[idx], XI = XI, BETA = BETA,
    cluster = TC %>% set_names(paste0("c", 1:m.)),
    cluster.labels = cluster, SNR = SNR,
    SDC = map(E, ~{
      k <- .x[1]; l <- .x[2]
      return(data.frame(
        pair = paste0(k, "-", l),
        score = sum((XI[k,] - XI[l,])^2)/sig2
      ))
    }) %>% exec(rbind, !!!.)
  ))
}
