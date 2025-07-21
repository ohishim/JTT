#' @title Consistent clustering for group-wise linear regression
#' @description \code{JTT} By join-two-together method, this function performs
#'     clustering for group-wise linear regression with graph structure (v0.1.0)
#'
#' @useDynLib JTT, .registration = TRUE
#'
#' @importFrom dplyr arrange pull
#' @importFrom magrittr %>% equals extract not set_names set_rownames
#' @importFrom purrr exec imap_dfr map map_chr map_dbl map_dfr map_lgl transpose
#' @importFrom stats optimize
#'
#' @param y a vector of a response variable
#' @param X a matrix of explanatory variables without intercept (intercept is
#'     automatically included)
#' @param group a vector of group indexes
#' @param adj a matrix with two columns expressing adjacency among groups
#' @param PSE post selection estimation; `"OLS"` (default) or `"penalized"`
#' @param alpha a penalty paramer of GCp criterion for JTT score
#' @param pse.alpha a penalty parameter of GCp criterion for post selection
#'     estimation when `PSE = "penalized"`
#' @param Rcpp if `TRUE`, C++ codes are used
#'
#' @return a list with the following elements:
#' \item{coefficients}{a matrix of estimated coefficients}
#'
#' \item{fitted.values}{a vector of fitted values}
#'
#' \item{cluster}{a list expressing group indexes of each cluster}
#'
#' \item{cluster.labels}{a vector of cluster indexes}
#'
#' \item{JTT.scores}{a dataframe with pairs of groups and JTT scores}
#'
#' \item{PostSelEst}{information of post selection estimation}
#'
#' \item{summary}{summary of estimation results}
#'
#' @export
#' @examples
#' #JTT(y, X, group, adj)

JTT <- function(y, X, group, adj, PSE="OLS", alpha=NULL, pse.alpha=NULL, Rcpp=FALSE){

  t1 <- proc.time()[3]

  if(Rcpp)
  {
    chol_solve <- chol_solve2
  } else
  {
    chol_solve <- chol_solve1
  }

  n <- length(y); m <- unique(group) %>% length; p <- ncol(X) + 1
  mp <- m*p; N <- n - mp

  X <- cbind(intercept=1, X)
  E <- apply(adj, 1, sort) %>% t %>% unique
  q <- nrow(E)

  id_ <- split(1:n, group)
  y_ <- split(y, group)
  X_ <- as.data.frame(X) %>% split(group) %>% map(data.matrix)

  res.mapm <- map(1:m, ~{
    Xj <- X_[[.x]]; yj <- y_[[.x]]
    M <- crossprod(Xj); X.y <- crossprod(Xj, yj) %>% drop
    y.Py <- crossprod(X.y, chol_solve(M, X.y)) %>% drop

    return(list(
      M = M, X.y = X.y, y.Py = y.Py,
      D0 = drop(crossprod(yj)) - y.Py
    ))
  }) %>% purrr::transpose()

  M <- res.mapm$M
  X.y <- res.mapm$X.y
  y.Py <- unlist(res.mapm$y.Py)
  D0 <- N/sum(unlist(res.mapm$D0))

  n0 <- map_dbl(y_, length) %>% min

  if(is.null(alpha))
  {
    B <- N*sqrt(N+p-2)/((N-2)*sqrt(N-4))
    beta <- B*(m^(1/4))*log(n0)/sqrt(p)
    alpha <- (N/(N - 2)) + beta
  }

  alp <- alpha*p

  if(Rcpp)
  {
    score <- JTT_score(q, m, E-1L, X.y, y.Py, M, D0, alp)
  } else
  {
    score <- map_dbl(1:q, function(j){
      Ej <- E[j,]; k <- Ej[1]; l <- Ej[2]
      X.y_kl <- X.y[[k]] + X.y[[l]]

      D1 <- sum(y.Py[Ej]) - drop(crossprod(
        X.y_kl, chol_solve(M[[k]]+M[[l]], X.y_kl)
      ))

      return(D0*D1 - alp)
    })
  }

  if(all(score > 0))
  {
    G <- as.list(1:m)
  } else
  {
    G <- create.cluster(E[score <= 0,,drop=F], m)
  }

  m. <- length(G)

  # cluster <- map(1:length(G), ~cbind(c=.x, g=G[[.x]])) %>% exec(rbind, !!!.) %>%
  #   as.data.frame %>% arrange(g) %>% use_series(c)
  cluster <- imap_dfr(G, ~data.frame(c=.y, g=.x)) %>% arrange(g) %>% pull(c)

  if(PSE == "OLS" | m. == 1)
  {
    Coef <- map(G, ~drop(
      chol_solve(M[.x] %>% Reduce(`+`, .), X.y[.x] %>% Reduce(`+`, .))
    )) %>% exec(rbind, !!!.) %>% magrittr::extract(cluster,)

    if(m. == 1)
    {
      PSE <- list(
        method = "OLS",
        warning = "only one cluster"
      )
    }
  } else if(PSE == "penalized")
  {
    PSEres <- penPSE(
      y, X,
      group = cluster[group],
      adj = cbind(cluster[adj[,1]], cluster[adj[,2]]) %>% unique %>%
        magrittr::extract(.[,1] - .[,2] != 0,),
      alpha = pse.alpha
    )

    Coef <- PSEres$coefficients %>% magrittr::extract(cluster,)

    PSE <- list(
      method = "penalized",
      lambda = PSEres$lambda,
      alpha = PSEres$alpha
    )
  }

  # Fit <- map(1:m, ~drop(X_[[.x]]%*%Coef[.x,])) %>% unlist %>%
  #   data.frame(fit=., id=unlist(id_)) %>% arrange(id) %>% dplyr::pull(fit)
  Fit <- map_dfr(1:m, ~data.frame(fit=drop(X_[[.x]]%*%Coef[.x,]), id=id_[[.x]])) %>%
    arrange(id) %>% pull(fit)

  t2 <- proc.time()[3]

  rss <- mean((y - Fit)^2)

  return(list(
    coefficients = Coef %>% set_rownames(paste0("c", cluster)),
    fitted.values = Fit,
    cluster = G %>% set_names(paste0("c", 1:m.)),
    cluster.labels = cluster %>% set_names(paste0("g", 1:m)),
    JTT.scores = data.frame(
      pair = map_chr(1:q, ~paste0(E[.x, 1], "-", E[.x, 2])),
      score = score
    ),
    PostSelEst = PSE,
    summary = data.frame(
      rss = rss,
      R2 = 1 - (rss / mean((y - mean(y))^2)),
      cluster = m.,
      alpha = alpha,
      runtime = (t2 - t1) %>% set_names(NULL)
    )
  ))
}
