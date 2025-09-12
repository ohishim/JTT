
chol_solve1 <- function(M, b){

  R <- chol(M)
  y <- backsolve(R, b, transpose=TRUE)
  return(backsolve(R, y))
}

corM <- function(k, rho){

  A <- map(1:k, ~{rho^abs(.x - (1:k))}) %>% exec(rbind, !!!.)

  Asq <- expm::sqrtm(A)

  return(list(
    M = A, sqrtM = Asq
  ))
}

create.cluster <- function(idx, m){

  G <- make_empty_graph(n = m, directed = FALSE) %>% add_edges(as.vector(t(idx)))

  labs <- components(G)$membership

  return(list(
    cluster = split(1:m, labs),
    labels = labs
  ))
}

true.cluster <- function(m, m.){

  if(m == 10)
  {
    if(m. == 3)
    {
      Tcluster <- list(c(1,2,3),c(4,5,6,9,10),c(7,8))
    }

    if(m. == 6)
    {
      Tcluster <- list(c(1,3),2,c(4,6,10),5,c(7,8),9)
    }
  } else if(m == 20)
  {
    if(m. == 6)
    {
      Tcluster <- list(c(1,2,3),c(4,5,6),c(7,8,19,20),c(9,10,12,13),c(11,14,15,16),c(17,18))
    }

    if(m. == 12)
    {
      Tcluster <- list(1,c(2,3),4,c(5,6),c(7,8),c(9,10),11,c(12,13),c(14,15,16),c(17,18),19,20)
    }
  } else if(m == 50)
  {
    if(m. == 15)
    {
      Tcluster <- list(
        c(1, 2, 3, 8), c(4, 14), c(5, 16, 17),
        c(6, 7, 18, 21, 22, 23), c(9, 10, 11, 13),
        c(12, 31, 40), c(15, 26, 27, 33, 49),  c(19, 20),
        c(24, 35, 45, 46), c(25, 30, 39, 44),
        c(28, 47, 48, 50), c(29, 37, 38), 32, c(34, 41), c(36, 42, 43)
      )
    }

    if(m. == 30)
    {
      Tcluster <- list(
        c(1, 2, 8), 3, c(4, 14), c(5, 16), c(6, 7), 9,
        c(10, 11), c(12, 31), 13, c(15, 26, 27, 33), c(17, 19), c(18, 21),
        20, c(22, 23), c(24, 45), c(25, 37), c(28, 47), c(29, 35, 38),
        c(30, 39), 32, 34, c(36, 43),
        40, 41, 42, 44, 46, 48, 49, 50
      )
    }
  }

  return(Tcluster %>% set_names(paste0("g", 1:m.)))
}

penPSE <- function(y, X, group, adj, chol_solve, alpha=NULL){

  m <- unique(group) %>% length

  y_ <- split(y, group)
  X_ <- as.data.frame(X) %>% split(group) %>% map(data.matrix)

  ns <- map_dbl(y_, length)

  LSE <- map(1:m, ~{
    chol_solve(
      crossprod(X_[[.x]]),
      crossprod(X_[[.x]], y_[[.x]]) %>% drop
    )
  }) %>% exec(rbind, !!!.)

  PSE <- map(1:m, ~{
    Dj <- adj[adj[,1] == .x, 2]
    Bj <- LSE[Dj,,drop=F]
    wj <- 1/sqrt(colSums((LSE[.x,] - t(Bj))^2))

    out <- .penPSE(
      y_[[.x]], X_[[.x]],
      colSums(wj*Bj)/sum(wj),
      alpha = alpha
    )
    return(out)
  }) %>% purrr::transpose()

  return(list(
    coefficients = PSE$coefficients %>% exec(rbind, !!!.),
    lambda = PSE$lambda %>% unlist,
    alpha = PSE$alpha %>% unlist
  ))
}

.penPSE <- function(y, X, b, alpha=NULL){

  n <- length(y); k <- ncol(X)

  if(is.null(alpha))
  {
    alpha <- 2*(n - k)/(n - k - 2)
  }

  svdX <- svd(X)
  d <- svdX$d^2
  z <- crossprod(svdX$u, y) %>% drop
  v <- crossprod(svdX$v, b) %>% drop %>% divide_by(svdX$d)
  s2 <- (sum(y^2) - sum(z^2))/(n-k)

  vec1 <- z - d*v

  f <- function(lam){
    map_dbl(lam, ~{
      del <- .x / (d + .x)
      return(
        (1/s2)*sum((del*vec1)^2) - alpha*sum(del)
      )
    })
  }

  lambda <- optimize(f, c(0, n*max(abs(y))))$minimum

  Coef <- drop(svdX$v %*% ((svdX$d/(d+lambda))*(z + lambda*v)))

  return(list(coefficients = Coef, lambda = lambda, alpha = alpha))
}

usethis::use_data(
  chol_solve1, corM, create.cluster, true.cluster, penPSE, .penPSE,
  internal=TRUE, overwrite=TRUE
)
