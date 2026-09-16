# helper functions for single-dataset meta analysis project
# Bayesian analysis
run_mama_bayes <- function(yi, sei, id_team, es_type, ni=NULL){

  if(es_type == "Cohen's d" | es_type == "Hedges' g"){
    outcome_scale <- sqrt(2)
  } else if(es_type == "log odds ratio" | es_type == "log hazard ratio" | es_type == "log relative risk"){
    outcome_scale <- sqrt(4)
  } else if(es_type == "z cor"){
    outcome_scale <- 1
  } else if(es_type == "beta"){
    if(is.null(ni)) stop("For standardized regression coefficients (beta), please provide the vector of sample sizes (ni) per team.")
    N <- median(ni)
    fit_scale     <- metafor::rma(yi = yi, sei = sei, method = "FE")
    outcome_scale <- fit_scale$se * sqrt(sum(length(yi) * N))
  } else {
    stop("The effect size type (es_type) needs to be one of: 'Cohen's d', 'Hedges' g', 'log odds ratio', 'log hazard ratio', 'log relative risk', 'z cor', or 'beta'.")
  }


  n_team <- length(unique(id_team))
  team_counts <- table(id_team)
  w <- mapply(function(team) {
    rep(1 / n_team / team_counts[team], team_counts[team])
  }, names(team_counts), SIMPLIFY = FALSE)
  w <- unlist(w)[order(match(id_team, names(team_counts)))]

  fit1 <- NoBMA(y = yi, se = sei,
                priors_effect        = prior("normal", list(mean = 0, sd = outcome_scale)),
                priors_heterogeneity = prior("normal", list(mean = 0, sd = outcome_scale/2), list(lower = 0, upper = Inf)),
                seed=2025)
  sfit1 <- summary(fit1, conditional = TRUE)
  tau_est_conditional <- sfit1$estimates_conditional["tau", "Median"]

  fit2 <- NoBMA(y = yi, se = sqrt(sei^2 / w),
                priors_effect             = prior("normal", list(mean = 0, sd = outcome_scale)),
                priors_heterogeneity      = prior("spike",  list(location = tau_est_conditional)),
                priors_heterogeneity_null = NULL, algorithm = "ss", adapt = 10000, burnin = 10000, sample = 10000,
                seed=2025)
  # I didn't implement bridge sampling for null model with fixed heterogeneity apparently (needs to switch to product space method)

  # summarize results
  return(list('normal' = fit1,
              'fractional' = fit2))
}

summarize_bayes <- function(mama){
  mu    <- median(mama$fractional$RoBMA$posteriors_conditional$mu)
  se    <- sd(mama$fractional$RoBMA$posteriors_conditional$mu)
  upper <- quantile(mama$fractional$RoBMA$posteriors_conditional$mu, p=.975)
  lower <- quantile(mama$fractional$RoBMA$posteriors_conditional$mu, p=.025)
  tau   <- median(mama$normal$RoBMA$posteriors_conditional$tau)
  tau_lower <- quantile(mama$normal$RoBMA$posteriors_conditional$tau, p=.025)
  tau_upper <- quantile(mama$normal$RoBMA$posteriors_conditional$tau, p=.975)
  evidence_mu  <- mama$fractional$RoBMA$inference['Effect','inclusion_BF']
  evidence_tau <- mama$normal$RoBMA$inference$Heterogeneity$BF
  return(data.frame(
    mu = mu, se = se, mu_upper = upper, mu_lower = lower,
    tau = tau, tau_lower = tau_lower, tau_upper = tau_upper,
    evidence_mu = evidence_mu, evidence_tau = evidence_tau
  ))
}


# frequentist analysis
run_mama <- function(yi, sei, id_team, es_type, ni=NULL){

  n_team <- length(unique(id_team))
  team_counts <- table(id_team)
  w <- mapply(function(team) {
    rep(1 / n_team / team_counts[team], team_counts[team])
  }, names(team_counts), SIMPLIFY = FALSE)
  w <- unlist(w)[order(match(id_team, names(team_counts)))]
  vi <- sei^2/w

  set.seed(2025)
  random1 <- try(metafor::rma(yi = yi, sei = sei, method = "REML"))
  random2 <- metafor::rma(yi = yi, vi = vi, tau2 = random1$tau2, method = "REML")

  confint_tau <- confint(random1)$random[2,2:3]
  bootstrapped_ci_tau <- FALSE

  if(confint_tau[2]==0){
    # bootstrap CI of tau
    B <- 2000
    tau_boot <- numeric(B)
    set.seed(2025)
    for(b in 1:B){
      # simulate study-level effects from fitted model
      # simulate sampling errors according to observed sei
      # predicted mean (mu_hat) can be used; here we use fitted values = overall intercept
      ystar <- rnorm(length(yi), mean = predict(random1)$pred, sd = sqrt(sei^2 + random1$tau2))
      fitb <- try(rma(yi = ystar, sei = sei, method="REML"), silent=TRUE)
      if(!inherits(fitb, "try-error")){
        tau_boot[b] <- fitb$tau2  # store tau^2 or take sqrt for tau
      } else {
        tau_boot[b] <- NA
      }
    }
    tau_boot <- na.omit(tau_boot)
    ci_tau2_boot <- quantile(tau_boot, probs = c(0.025, 0.975))
    ci_tau_boot <- sqrt(ci_tau2_boot)
    confint_tau <- ci_tau_boot
    bootstrapped_ci_tau = TRUE
  }

  if(es_type == "Cohen's d" | es_type == "Hedges' g"){
    outcome_scale <- sqrt(2)
  } else if(es_type == "log odds ratio" | es_type == "log hazard ratio" | es_type == "log relative risk"){
    outcome_scale <- sqrt(4)
  } else if(es_type == "z cor"){
    outcome_scale <- 1
  } else if(es_type == "beta"){
    if(is.null(ni)) stop("For standardized regression coefficients (beta), please provide the vector of sample sizes (ni) per team.")
    N <- median(ni)
    fit_scale     <- metafor::rma(yi = yi, sei = sei, method = "FE")
    outcome_scale <- fit_scale$se * sqrt(sum(length(yi) * N))
  } else {
    stop("The effect size type (es_type) needs to be one of: 'Cohen's d', 'Hedges' g', 'log odds ratio', 'log hazard ratio', 'log relative risk', 'z cor', or 'beta'.")
  }

  return(list('normal' = data.frame(
    mu = random1$b,
    se = random1$se,
    upper = random1$ci.ub,
    lower = random1$ci.lb,
    pval = random1$pval,
    tau2 = random1$tau2,
    Qtau2 = random1$QE,
    pval_tau2 = random1$QEp,
    ci_tau_lower = confint_tau[1],
    ci_tau_upper = confint_tau[2],
    bootstrapped_ci_tau = bootstrapped_ci_tau
  ),
  'fractional' = data.frame(
    mu = random2$b,
    se = random2$se,
    upper = random2$ci.ub,
    lower = random2$ci.lb,
    pval = random2$pval,
    tau2 = random2$tau2
  ),
  n_team = n_team,
  i2 = random1$I2
  ))
}

summarize_metafor <- function(mama){
  mu <- mama$fractional$mu
  se <- mama$fractional$se
  upper <- mama$fractional$upper
  lower <- mama$fractional$lower
  tau <- sqrt(mama$normal$tau2)
  evidence_mu <- mama$fractional$bf10
  evidence_tau <- ifelse(mama$normal$pval_tau2 < .05, "significant", "not significant")
  return(data.frame(
    mu = mu, se = se, upper = upper, lower = lower,
    tau = tau, evidence_mu = evidence_mu, evidence_tau = evidence_tau
  ))
}


# plotting functions
forest_mama <- function(x1, x2, slab = "Study", es_type = "beta", transform=FALSE,
                        title = NULL, subtitle = NULL, xmax = NULL, large=FALSE,
                        labels = NULL, digits = 2, x_labels = NULL, xmin = NULL,
                        est_lab = NULL){
  x2a <- x2$normal
  sum <- summarize_bayes(x2)
  est_mu <- c(x1$normal$mu, median(x2a$RoBMA$posteriors_conditional$mu), x1$fractional$mu, sum$mu)
  lCI_mu <- c(x1$normal$lower, quantile(x2a$RoBMA$posteriors_conditional$mu, p=0.025), x1$fractional$lower, sum$mu_lower)
  uCI_mu <- c(x1$normal$upper, quantile(x2a$RoBMA$posteriors_conditional$mu, p=0.975), x1$fractional$upper, sum$mu_upper)
  data=x2a[["data"]]
  data <- data.frame(data)
  if (length(slab) == 1) {
    data$study <- paste0(slab, " ", seq_len(nrow(data)))
  } else if (length(slab) == nrow(data)) {
    data$study <- slab
  } else {
    stop(
      "`slab` must be either one label or a vector with one label per estimate."
    )
  }
  data$lower <- data$y - 1.96 * data$se
  data$upper <- data$y + 1.96 * data$se
  if(transform) {
    data_exp <- data
    data_exp$y <- exp(data_exp$y)
    data_exp$lower <- exp(data_exp$lower)
    data_exp$upper <- exp(data_exp$upper)
    est_mu_exp <- exp(est_mu)
    lCI_mu_exp <- exp(lCI_mu)
    uCI_mu_exp <- exp(uCI_mu)
    x_breaks <- log(x_labels)
    data_exp <- data_exp %>% arrange(desc(y))
  }
  data <- data[order(data$y, decreasing = TRUE), ]
  data$x <- (nrow(data) + 5):6

  idx <- c(1, nrow(data))
  y_at <- c(1:4, data$x)
  if(large) y_at <- c(c(-0.17, -.13, -.09, -.05) * nrow(data), data$x[idx])
  y_labels <- c("Unadjusted Classical", "Unadjusted Bayesian",
                "Adjusted Classical", "Adjusted Bayesian",
                data$study)
  y_labels2 <- paste0(format(round(c(est_mu, data$y), digits), nsmall = 2),
                      " [", format(round(c(lCI_mu, data$lower), digits), nsmall = 2),
                      ", ", format(round(c(uCI_mu, data$upper), digits), nsmall = 2),
                      "]")
  if(transform) y_labels2 <- paste0(format(round(c(est_mu_exp, data_exp$y), digits), nsmall = 2),
                                    " [", format(round(c(lCI_mu_exp, data_exp$lower), digits), nsmall = 2),
                                    ", ", format(round(c(uCI_mu_exp, data_exp$upper), digits), nsmall = 2),
                                    "]")
  if(large) y_labels <- c("Unadjusted Classical", "Unadjusted Bayesian",
                          "Adjusted Classical", "Adjusted Bayesian", data$study[idx])
  if(large) y_labels2 <- c(paste0(format(round(c(est_mu), digits), nsmall = 2),
                                  " [", format(round(c(lCI_mu), digits), nsmall = 2),
                                  ", ", format(round(c(uCI_mu), digits), nsmall = 2),
                                  "]"), "", est_lab)
  if(large&transform) y_labels2 <- c(paste0(format(round(c(est_mu_exp), digits), nsmall = 2),
                                            " [", format(round(c(lCI_mu_exp), digits), nsmall = 2),
                                            ", ", format(round(c(uCI_mu_exp), digits), nsmall = 2),
                                            "]"), "", est_lab)
  ylim <- c((y_at[1]-1), max(data$x) + 1)
  if(large) ylim <- c((y_at[1]-0.008 * nrow(data)), max(data$x) + 1)
  xlab <- es_type
  if(es_type=="beta") xlab <- expression("Standardized " * beta)
  if(!transform){
    x_breaks <- pretty(range(c(data$lower, data$upper, lCI_mu, uCI_mu)))
    x_labels <- x_breaks
  }
  xlim <- range(x_breaks)

  #if(!is.null(xmax)){
  #  data$upper[data$upper > xmax] <- xmax
  #  data$flag_upper <- ifelse(data$upper == xmax, 1, 0)
  #}
  #if(!is.null(xmin)){
  #  data$lower[data$lower < xmin] <- xmin
  #  data$flag_lower <- ifelse(data$lower == xmin, 1, 0)
  #}

  x_labels <- c(paste0(x_labels[1],"\n",labels[1]),  x_labels[2:(length(x_labels) - 1)],
                              paste0(x_labels[length(x_labels)], "\n",labels[2]))
  alpha <- ifelse(large, 0.1, 1)
  size <- ifelse(large, 0.8, 1)
  y_comb <- c(1, 1.25, 1, 0.75, 2, 2.25, 2, 1.75, 3, 3.25, 3, 2.75, 4, 4.25, 4, 3.75)
  #if(large) y_comb <- c(y_at[1], y_at[1]+0.007 * nrow(data), y_at[1], y_at[1]-0.007 * nrow(data))
  if(large) y_comb <- as.vector(sapply(y_at[1:4], function(y) c(y, y + 0.007 * nrow(data), y, y - 0.007 * nrow(data))))

  # Split data into regular and truncated
  data$flag_upper <- FALSE
  data$flag_lower <- FALSE

  # Truncate and flag upper bounds if xmax specified
  if (!is.null(xmax)) {
    data$flag_upper <- data$upper > xmax
    data$upper[data$flag_upper] <- xmax
  }

  # Truncate and flag lower bounds if xmin specified
  if (!is.null(xmin)) {
    data$flag_lower <- data$lower < xmin
    data$lower[data$flag_lower] <- xmin
  }

  # Split datasets
  data_regular <- subset(data, !flag_upper & !flag_lower)
  data_right_only <- subset(data, flag_upper & !flag_lower)
  data_left_only  <- subset(data, flag_lower & !flag_upper)
  data_both_sides <- subset(data, flag_upper & flag_lower)

  df_polys <- data.frame(
    x = as.vector(rbind(lCI_mu, est_mu, uCI_mu, est_mu)),
    y = y_comb,
    group = rep(1:4, each = 4)
  )

  plot <- ggplot2::ggplot()
  #plot <- plot + ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = data$lower,
  #                                                              xmax = data$upper, y = data$x),
  #                                       color = "black", alpha = alpha, height = 0.25)
  # Normal errorbars
  plot <- plot +
    ggplot2::geom_errorbarh(
      data = data_regular,
      mapping = ggplot2::aes(xmin = lower, xmax = upper, y = x),
      color = "black", alpha = alpha, height = 0.25
    )

  # Add right-arrow for truncated upper bounds
  if (nrow(data_right_only) > 0) {
    plot <- plot +
      ggplot2::geom_segment(
        data = data_right_only,
        mapping = ggplot2::aes(x = lower, xend = xmax, y = x, yend = x),
        arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "open"),
        color = "black", alpha = alpha
      ) +
      ggplot2::geom_segment(
        data = data_right_only,
        mapping = ggplot2::aes(x = lower, xend = lower, y = x - 0.1, yend = x + 0.1),
        color = "black", alpha = alpha
      )
  }

  # Add left-arrow for truncated lower bounds
  if (nrow(data_left_only) > 0) {
    plot <- plot +
      ggplot2::geom_segment(
        data = data_left_only,
        mapping = ggplot2::aes(x = xmin, xend = upper, y = x, yend = x),
        arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "open", ends = "first"),
        color = "black", alpha = alpha
      ) +
      ggplot2::geom_segment(
        data = data_left_only,
        mapping = ggplot2::aes(x = upper, xend = upper, y = x - 0.1, yend = x + 0.1),
        color = "black", alpha = alpha
      )
  }
  if (nrow(data_both_sides) > 0) {
    plot <- plot +
      ggplot2::geom_segment(
        data = data_both_sides,
        mapping = ggplot2::aes(x = lower, xend = upper, y = x, yend = x),
        arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "open", ends = "both"),
        color = "black", alpha = alpha
      )
  }
  plot <- plot + ggplot2::geom_point(mapping = ggplot2::aes(x = data$y,
                                                            y = data$x), shape = 15, size = size)
  plot <- plot + ggplot2::geom_polygon(data= df_polys,
                                       mapping = ggplot2::aes(x = x,y = y, group=group),
                                       fill = "black")
  plot <- plot + ggplot2::geom_line(mapping = ggplot2::aes(x = c(0,
                                                                 0), y = ylim), linetype = "dotted")

  plot <- plot + ggplot2::scale_y_continuous(name = "",
                                             breaks = y_at, labels = y_labels, limits = ylim,
                                             sec.axis = ggplot2::sec_axis(~., breaks = y_at, labels = y_labels2))
  plot <- plot + ggplot2::scale_x_continuous(name = xlab,
                                             breaks = x_breaks, labels = x_labels, limits = xlim)
  plot <- plot + ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                                axis.line.y = ggplot2::element_blank(),
                                axis.ticks.y = ggplot2::element_blank(),
                                axis.line.x = ggplot2::element_line(arrow = grid::arrow(length = unit(0.2, "cm"),
                                                                                        ends = "both")),
                                axis.ticks.x = ggplot2::element_line(linewidth = 0.5),
                                axis.text.y = ggplot2::element_text(hjust = 0, color = "black"),
                                axis.text.y.right = ggtext::element_markdown(hjust = 0.5,
                                                                          color = "black"))
  if(large) plot <- plot + base_breaks_y(c(1,nrow(data)))
  #plot <- plot + ggplot2::ggtitle(title, subtitle = subtitle)
  p <- plot
  if(!large){
    p <- ggdraw(plot) +
      draw_label(est_lab,
                 x = 0.96, y = 0.96,
                 angle = 0, hjust = 1, vjust = 1,
                 size = 9, fontface = "bold")
  }
  return(p)
}

# define custom y axis function
base_breaks_y <- function(yLimits) {
  d <- data.frame(x=-Inf, xend=-Inf, y=yLimits[1], yend=yLimits[2])
  list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y, xend=xend, yend=yend), size = 0.75, inherit.aes=FALSE))
}
base_breaks_x <- function(xLimits) {
  d <- data.frame(x=xLimits[1], xend=xLimits[2], y=-Inf, yend=-Inf)
  list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y, xend=xend, yend=yend), size = 0.75, inherit.aes=FALSE))
}


# function to extract key values from one SDMA for summary
extract_sdma_summary <- function(
    sdma_name,
    predictor,
    estimand,
    contrast,
    data,
    mama,
    mama_bayes
) {
  
  # Fractionally adjusted Bayesian summary
  bayes_summary <- summarize_bayes(mama_bayes)
  
  # Unadjusted conditional Bayesian posterior for mu
  unadjusted_mu <- quantile(
    mama_bayes$normal$RoBMA$posteriors_conditional$mu,
    probs = c(0.025, 0.50, 0.975),
    na.rm = TRUE
  )
  
  # Adjusted Bayesian inclusion BF for the pooled effect
  bayes_adjusted_BF10_effect <-
    mama_bayes$fractional$RoBMA$inference[
      "Effect",
      "inclusion_BF"
    ]
  
  # Unadjusted Bayesian inclusion BF for the pooled effect
  bayes_unadjusted_BF10_effect <-
    mama_bayes$normal$RoBMA$inference$Effect$BF
  
  tibble(
    sdma = sdma_name,
    predictor = predictor,
    estimand = estimand,
    contrast = contrast,
    
    # Number of included estimates and teams
    n_estimates = nrow(data),
    n_teams = n_distinct(data$team),
    
    # Distribution of analysis-specific estimates
    median_estimate = median(data$y, na.rm = TRUE),
    minimum_estimate = min(data$y, na.rm = TRUE),
    maximum_estimate = max(data$y, na.rm = TRUE),
    
    # Fractionally adjusted frequentist results
    freq_adjusted_mu = as.numeric(mama$fractional$mu),
    freq_adjusted_ci_lower = as.numeric(mama$fractional$lower),
    freq_adjusted_ci_upper = as.numeric(mama$fractional$upper),
    freq_adjusted_p = as.numeric(mama$fractional$pval),
    
    # Unadjusted frequentist results
    freq_unadjusted_mu = as.numeric(mama$normal$mu),
    freq_unadjusted_ci_lower = as.numeric(mama$normal$lower),
    freq_unadjusted_ci_upper = as.numeric(mama$normal$upper),
    freq_unadjusted_p = as.numeric(mama$normal$pval),
    
    # Frequentist heterogeneity
    freq_tau = sqrt(as.numeric(mama$normal$tau2)),
    freq_heterogeneity_p = as.numeric(mama$normal$pval_tau2),
    
    # Fractionally adjusted Bayesian pooled effect
    bayes_adjusted_mu = as.numeric(bayes_summary$mu),
    bayes_adjusted_ci_lower = as.numeric(bayes_summary$mu_lower),
    bayes_adjusted_ci_upper = as.numeric(bayes_summary$mu_upper),
    
    # Fractionally adjusted Bayesian evidence for the effect
    bayes_adjusted_BF10_effect =
      as.numeric(bayes_adjusted_BF10_effect),
    
    bayes_adjusted_BF01_effect =
      as.numeric(1 / bayes_adjusted_BF10_effect),
    
    # Fractionally adjusted Bayesian heterogeneity
    bayes_tau = as.numeric(bayes_summary$tau),
    bayes_tau_ci_lower = as.numeric(bayes_summary$tau_lower),
    bayes_tau_ci_upper = as.numeric(bayes_summary$tau_upper),
    bayes_BF10_heterogeneity =
      as.numeric(bayes_summary$evidence_tau),
    
    bayes_BF01_heterogeneity =
      as.numeric(1 / bayes_summary$evidence_tau),
    
    # Unadjusted conditional Bayesian pooled effect
    bayes_unadjusted_mu_median = as.numeric(unadjusted_mu["50%"]),
    bayes_unadjusted_ci_lower = as.numeric(unadjusted_mu["2.5%"]),
    bayes_unadjusted_ci_upper = as.numeric(unadjusted_mu["97.5%"]),
    
    # Unadjusted Bayesian evidence for the effect
    bayes_unadjusted_BF10_effect =
      as.numeric(bayes_unadjusted_BF10_effect),
    
    bayes_unadjusted_BF01_effect =
      as.numeric(1 / bayes_unadjusted_BF10_effect)
  )
}
