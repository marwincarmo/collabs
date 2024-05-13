
# --- Expected Influence --- 

get_ei <- function(net, data = NULL, type = "pcor") {
  # net <- BGGM::estimate(Y, progress = FALSE)
  
  if (type == "pcor") {
    post_samps <- BGGM::posterior_samples(net)
    cnames <- colnames(net$Y)
  } else if (type == "cor") {
    post_samps <- posterior_samples_zero_order(net, data = data)
    cnames <- colnames(data)
  } else {
    stop("Check type")
  }
  
  cnamesPost <- colnames(post_samps)
  
  post_ei <- array(NA,
                   dim = c(nrow(post_samps), length(cnames)),
                   dimnames = list(NULL, c(cnames)))
  
  calc_ei <- function(node) {
    edge_idx <- grep(node, cnamesPost)
    current_nodes <- post_samps[, edge_idx]
    ei <- rowSums(current_nodes)
    post_ei[, node] <- ei
  }
  
  post <- sapply(cnames, 
                 function(x) calc_ei(x))
  return(post)
}


get_ei_diffs <- function(net1, net2, cred = 0.95, data = NULL, type = "pcor") {
  lb <- (1 - cred) / 2
  ub <- 1 - lb
  
  # summary group 1
  ei1 <- get_ei(net1, data = data, type = type)
  m1 <- colMeans(ei1)
  sd_1 <- apply(ei1, 2, sd)
  lower_1 <- apply(ei1, 2, function(x) quantile(x, lb))
  upper_1 <- apply(ei1, 2, function(x) quantile(x, ub))
  df1 <- round(data.frame(Post.mean = m1, SD = sd_1, cred.lb = lower_1, cred.ub = upper_1), 3)
  # summary group2
  ei2 <- get_ei(net2, data = data, type = type)
  m2 <- colMeans(ei2)
  sd_2 <- apply(ei2, 2, sd)
  lower_2 <- apply(ei2, 2, function(x) quantile(x, lb))
  upper_2 <- apply(ei2, 2, function(x) quantile(x, ub))
  df2 <- round(data.frame(Post.mean = m2, SD = sd_2, cred.lb = lower_2, cred.ub = upper_2), 3)
  # summary differences
  diffs <- ei1 - ei2
  
  means <- colMeans(diffs)
  sd_diff <- apply(diffs, 2, sd)
  lower_q <- apply(diffs, 2, function(x) quantile(x, lb))
  upper_q <- apply(diffs, 2, function(x) quantile(x, ub))
  
  totals <- round(data.frame(Post.mean = means, SD = sd_diff, cred.lb = lower_q, cred.ub = upper_q), 3)
  
  all <- list(Yg1 = df1, Yg2 = df2, `EI Diff (Yg1 - Yg2)` = totals)
  return(all)
}



# --- Global Strength ---

get_gs <- function(net, cred = 0.95, summary = FALSE, type = "pcor", data = NULL) {
  
  lq <- (1 - cred) / 2
  uq <- 1 - lq
  
  if (type == "pcor") {
    post_samps <- BGGM::posterior_samples(net)
  } else if (type == "cor") {
    post_samps <- posterior_samples_zero_order(net, data = data)
  } else {
    stop("Check type")
  }
  
  lb <- apply(post_samps, 2, function(x) quantile(x, lq))
  ub <- apply(post_samps, 2, function(x) quantile(x, uq))
  
  
  excludes_zero <- lb >= 0  | ub <= 0
  post_samps_zero <- t(apply(post_samps, 1, function(x) x * excludes_zero))
  gs <- rowSums(abs(post_samps_zero))
  
  return(gs)
}

# data arg not needed for gs, but for posterior samples for zero order cor
get_gs_diffs <- function(net1, net2, cred = 0.95, type = "pcor", data = NULL) {
  lb <- (1 - cred) / 2
  ub <- 1 - lb
  
  # summary for group 1
  gs1 <- get_gs(net1, type = type, data = data)
  m1 <- mean(gs1)
  sd_1 <- sd(gs1)
  lower_1 <-  quantile(gs1, lb)
  upper_1 <- quantile(gs1, ub)
  
  # summary for group 2
  gs2 <- get_gs(net2, type = type, data = data)
  m2 <- mean(gs2)
  sd_2 <- sd(gs2)
  lower_2 <-  quantile(gs2, lb)
  upper_2 <- quantile(gs2, ub)
  
  # summary for differences
  diffs <- abs(gs1 - gs2)
  
  mean_diff <- mean(diffs)
  sd_diff <- sd(diffs)
  lower_diff <-  quantile(diffs, lb)
  upper_diff <- quantile(diffs, ub)
  
  # concatenate summaries
  means <- c(m1, m2, mean_diff)
  sds <- c(sd_1, sd_2, sd_diff)
  ci.lb <- c(lower_1, lower_2, lower_diff)
  ci.ub <- c(upper_1, upper_2, upper_diff)
  
  totals <- round(data.frame(Post.mean = means, SD = sds, cred.lb = ci.lb, cred.ub = ci.ub), 3)
  rownames(totals) <- c("Group1", "Group2", "GS Diff")
  return(totals)
}

f_gs_ppc <- function(x, y) {
  
  fit1 <-  estimate(x, progress = F)
  fit2 <-  estimate(y, progress = F)
  
  gs_diffs <- get_gs_diffs(fit1, fit2, type = type, data = data)
  
  
  return(gs_diffs["GS Diff", "Post.mean"])
}

f_gs_cor <- function(x, y) {
  
  fit1 <-  zero_order_cors(x, progress = F)
  fit2 <-  zero_order_cors(y, progress = F)

  gs_diffs <- get_gs_diffs(fit1, fit2, type = "cor", data = NULL)
  
  return(gs_diffs["GS Diff", "Post.mean"])
}

# -- Make groups for age analyses --
make_groups <- function(data, grp, vars) {
  current_age <- 
    data %>% 
    dplyr::filter(age_grp == grp) %>% 
    dplyr::select(dplyr::all_of(vars))
  prior_age1 <- 
    data %>% 
    dplyr::filter(age_grp == (grp - 2)) %>% 
    dplyr::select(dplyr::all_of(vars))
  prior_age2 <- 
    data %>% 
    dplyr::filter(age_grp == (grp - 1)) %>% 
    dplyr::select(dplyr::all_of(vars))
  
  return(list(current_age = current_age, 
              prior_age1 = prior_age1,
              prior_age2 = prior_age2))
}

# -- Make groups for moderator analyses --
make_groups_adhd <- function(data, age_bin, vars, mod) {
  grouped <- data %>% 
    filter(age_grp == age_bin) %>% 
    group_by_(all_of(mod)) %>% 
    select(all_of(vars))
  
  split_groups <- group_split(grouped)
  split_groups_df <- lapply(split_groups, as.data.frame)
  split_groups_df <- lapply(split_groups_df, function(x) select(x, -all_of(mod)))
  return(split_groups_df)
}

# -- Get Descriptives --

get_descriptives <- function(net) {
  
  imputed <- sweep(net$Y, 2, attr(net$Y, "scaled:center"), "+")
  
  means <- colMeans(imputed)
  sds <- apply(imputed, 2, sd)
  
  lbs <- apply(imputed, 2, function(x) quantile(x, 0.025))
  ubs <- apply(imputed, 2, function(x) quantile(x, 0.975))
  
  desc <- data.frame(Mean = means, SD = sds, cred.lb = lbs, cred.ub = ubs)
  
  return(round(desc, 3))
}



#####################################
# -- Functions for zero order analyses --
#####################################
posterior_samples_zero_order <- function(object, data){
  R <- object$R
  
  Rdims <- dim(R)
  # nodes
  p <- Rdims[[2]]
  
  # total partials
  cors_total <- p * (p - 1) * 0.5
  
  # identity matrix
  I_p <- diag(p)
  
  # iterations
  iter <- Rdims[[3]]
  
  # pcor samples
  pcor_samples <-
    matrix(
      R[, , 1:iter][upper.tri(I_p)],
      nrow =  iter,
      ncol = cors_total,
      byrow = TRUE
    )
  
  # column names
  cn <- colnames(data)
  
  if(is.null(cn)){
    
    col_names <- sapply(1:p, function(x)  paste(1:p, x, sep = "--"))[upper.tri(I_p)]
    
  } else {
    
    col_names <- sapply(cn, function(x)  paste(cn, x, sep = "--"))[upper.tri(I_p)]
  }
  
  
  colnames(pcor_samples) <- col_names
  
  posterior_samples <- pcor_samples
  
  
  return(posterior_samples)
  
}

# -- get edges for plotting in zero order --

get_zero_order_diffs <- function(net1, net2, cred = 0.95, data = data, summary = FALSE) {
 
 lb <- (1 - cred)/2
 ub <- 1 - lb
 
 post1 <- posterior_samples_zero_order(net1, data = data) 
 post2 <- posterior_samples_zero_order(net2, data = data) 
 
 # summary for group 1
 m1 <- colMeans(post1)
 sd_1 <- apply(post1, 2, sd)
 lwr_1 <-  apply(post1, 2, quantile, lb)
 upr_1 <- apply(post1, 2, quantile, ub)
 
 # summary for group 2
 m2 <- colMeans(post2)
 sd_2 <- apply(post2, 2, sd)
 lwr_2 <-  apply(post2, 2, quantile, lb)
 upr_2 <- apply(post2, 2, quantile, ub)
 
 # summary for diff
 diff <- post1 - post2
 mdiff <- colMeans(diff)
 sd_diff <- apply(diff, 2, sd)
 lwr_diff <-  apply(diff, 2, quantile, lb)
 upr_diff <- apply(diff, 2, quantile, ub)
 
 if (summary) {
   d1 <- data.frame(Post.mean = m1, SD = sd_1, cred.lb = lwr_1, cred.ub = upr_1)
   d2 <- data.frame(Post.mean = m2, SD = sd_2, cred.lb = lwr_2, cred.ub = upr_2)
   
   return(list(d1, d2))
 }
 
 return(data.frame(Mean = mdiff, SD =  sd_diff, cred.lb = lwr_diff, cred.ub = upr_diff))
 
}

# -- plot differences in edges --

plot_edges_diff <- function(net1, net2, cred = 0.95, data = data) {
  
  diff_data <- get_zero_order_diffs(net1, net2, cred, data)
  diff_data_ordered <- diff_data[order(diff_data$Mean, decreasing = FALSE), ]
  diff_data_ordered$Relation <- factor(rownames(diff_data_ordered),
                                       levels = rownames(diff_data_ordered),
                                       labels = rownames(diff_data_ordered))
  ggplot(diff_data_ordered, aes(x = Relation, y = Mean)) +
    geom_point(size = 2) +
    xlab("Index") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    geom_errorbar(aes(ymax = cred.ub, ymin = cred.lb), width = 0, color = "black") +
    ggtitle("Yg1 - Yg2")  -> p
  
  return(p)
}
