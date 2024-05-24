roll_your_own <- function(object,
                          FUN,
                          iter = NULL,
                          select = FALSE,
                          cred = 0.95,
                          progress = TRUE,
                          ...){
  
  # if(! all( c("estimate", "default") %in% class(object)) ){
  #   stop("class must be 'estimate'")
  # }
  
  if(!isFALSE(select)){
    
    sel <- BGGM::select(object, cred = cred)
    adj <- sel$adj
    
  } else {
    
    p <- ncol(object$pcor_mat)
    adj <- matrix(1, p, p)
    
  }
  
  if(is.null(iter)){
    
    iter <- object$iter
    
  }
  
  pcors <- object$fit$pcors[, , 51:(iter + 50)]
  
  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }
  
  results <- sapply(1:iter, function(x) {
    
    pcors_s <- pcors[, , x] * adj
    
    est <- FUN(pcors_s, communities=1:21)
    
    if(isTRUE(progress)){
      
      utils::setTxtProgressBar(pb, x)
      
    }
    
    est
    
  })
  
  returned_object <- list(results = results, iter = iter)
  
  class(returned_object) <- c("BGGM",
                              "roll_your_own")
  
  return(returned_object)
}


print_roll_your_own <- function(x, cred = 0.95, ...) {
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Network Stats: Roll Your Own\n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("--- \n")
  cat("Estimates: \n\n")
  
  lb <- (1-cred) / 2
  ub <- 1 - lb
  
  dims <- dim(x$results)
  
  if(is.null(dims)){
    
    mu <- mean(x$results)
    
    scale <- sd(x$results)
    
    res <- data.frame(Post.mean = round(mean(x$results), 3),
                      Post.sd =    round(sd(x$results), 3),
                      Cred.lb = round(quantile(x$results, probs = lb), 3),
                      Cred.ub = round(quantile(x$results, probs = ub), 3) )
  } else {
    
    mu <-  apply( x$results, 1, mean)
    p <- length(mu)
    scale <- apply( x$results, 1, sd)
    ci_lb <- apply( x$results, 1, quantile, lb)
    ci_ub <- apply( x$results, 1, quantile, ub)
    
    res<- data.frame(Node = 1:p,
                     Post.mean = round(mu, 3),
                     Post.sd = round(scale, 3),
                     Cred.lb = round(ci_lb, 3),
                     Cred.ub = round(ci_ub, 3))
  }
  
  print(res, row.names = FALSE)
  cat("--- \n")
}