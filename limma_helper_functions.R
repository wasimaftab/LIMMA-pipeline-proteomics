##########Function Blocks##################################################################################
eb.fit <- function(dat, design, gene) {
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma) ^ 2
  s2.post <- fit.eb$s2.post
  t.ord <-
    fit.eb$coefficients[, 2] / fit.eb$sigma / fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2 * pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <-
    data.frame(logFC,
               t.ord,
               t.mod,
               p.ord,
               p.mod,
               q.ord,
               q.mod,
               df.r,
               df.0,
               s2.0,
               s2,
               s2.post,
               gene)
  return(results.eb)
}

readinteger <- function(str)
{
  n <- readline(prompt = str)
  n <- as.integer(n)
  if (is.na(n) || (n <= 0)) {
    print("Enter positive integer only")
    n <- readinteger(str)
  }
  return(n)
}

readfloat <- function(str)
{
  n <- readline(prompt = str)
  n <- as.double(n)
  if (is.na(n) || (n <= 0)) {
    print("Enter a positive number only")
    n <- readfloat(str)
  }
  return(n)
}

## fitting a normal distribution
fitnormal <- function (x, exact = TRUE) {
  if (exact) {
    ################################################
    ## Exact inference based on likelihood theory ##
    ################################################
    ## minimum negative log-likelihood (maximum log-likelihood) estimator of `mu` and `phi = sigma ^ 2`
    n <- length(x)
    mu <- sum(x) / n
    phi <- crossprod(x - mu)[1L] / n  # (a bised estimator, though)
    ## inverse of Fisher information matrix evaluated at MLE
    invI <- matrix(c(phi, 0, 0, phi * phi), 2L,
                   dimnames = list(c("mu", "sigma2"), c("mu", "sigma2")))
    ## log-likelihood at MLE
    loglik <- -(n / 2) * (log(2 * pi * phi) + 1)
    ## return
    return(list(
      theta = c(mu = mu, sigma2 = phi),
      vcov = invI,
      loglik = loglik,
      n = n
    ))
  }
  else {
    ##################################################################
    ## Numerical optimization by minimizing negative log-likelihood ##
    ##################################################################
    ## negative log-likelihood function
    ## define `theta = c(mu, phi)` in order to use `optim`
    nllik <- function (theta, x) {
      (length(x) / 2) * log(2 * pi * theta[2]) + crossprod(x - theta[1])[1] / (2 * theta[2])
    }
    ## gradient function (remember to flip the sign when using partial derivative result of log-likelihood)
    ## define `theta = c(mu, phi)` in order to use `optim`
    gradient <- function (theta, x) {
      pl2pmu <- -sum(x - theta[1]) / theta[2]
      pl2pphi <-
        -crossprod(x - theta[1])[1] / (2 * theta[2] ^ 2) + length(x) / (2 * theta[2])
      c(pl2pmu, pl2pphi)
    }
    ## ask `optim` to return Hessian matrix by `hessian = TRUE`
    ## use "..." part to pass `x` as additional / further argument to "fn" and "gn"
    ## note, we want `phi` as positive so box constraint is used, with "L-BFGS-B" method chosen
    init <-
      c(sample(x, 1), sample(abs(x) + 0.1, 1))  ## arbitrary valid starting values
    z <-
      optim(
        par = init,
        fn = nllik,
        gr = gradient,
        x = x,
        lower = c(-Inf, 0),
        method = "L-BFGS-B",
        hessian = TRUE
      )
    ## post processing ##
    theta <- z$par
    loglik <- -z$value  ## flip the sign to get log-likelihood
    n <- length(x)
    ## Fisher information matrix (don't flip the sign as this is the Hessian for negative log-likelihood)
    I <- z$hessian / n  ## remember to take average to get mean
    invI <- solve(I, diag(2L))  ## numerical inverse
    dimnames(invI) <- list(c("mu", "sigma2"), c("mu", "sigma2"))
    ## return
    return(list(
      theta = theta,
      vcov = invI,
      loglik = loglik,
      n = n
    ))
  }
}


readinteger_binary <- function(str)
{
  n <- readline(prompt = str)
  n <- as.integer(n)
  # if (is.na(n) || (n <= 0)){
  if (is.na(n) || (n < 0) || (n > 1)) {
    print("Enter positive integer only")
    n <- readinteger_binary(str)
  }
  return(n)
}

data_sanity_check <- function(temp, exprement, exp_str) {
  exprement_reps <-
    # select(temp, matches(paste('.*', exp_str, '.*', sep = ''))) # regex to match substr ^.*[D][E][F].*$
    select(temp, matches(paste('^.*', exp_str, '.*$', sep = '')))
  row <- nrow(exprement_reps)
  col <- ncol(exprement_reps)
  if (!row * col || suppressWarnings(!is.na(as.integer(exp_str)))) {
    cat('Check', exprement, 'name and enter correct one\n')
    exp_str <-
      readline(
        cat(
          'Enter',
          exprement,
          'name(case insensitive) as it appeared in the iBAQ/LFQ column= '
        )
      )
    exprement_reps <- data_sanity_check(temp, 'treatment', exp_str)
  }
  return(exprement_reps)
}

display_plotly_figs <- function(dat, FC_Cutoff, filename_mod, filename_ord) {
  f <- list(family = "Arial, sans-serif",
            size = 18,
            color = "#7f7f7f")
  x_axis <- list(
    # title = "Log2 Fold Changes",
    constraintoward = "bottom",
    title = "log2 fold change",
    titlefont = f,
    showgrid = FALSE,
    showticklabels = TRUE,
    autotick = FALSE,
    ticks = "outside",
    tick0 = 0,
    dtick = 1,
    tickwidth = 1,
    position = -1,
    rangemode = "tozero"
  )
  y_axis <- list(
    # title = "-Log10 P-Value",
    title = "-log10 p-value",
    titlefont = f,
    showgrid = FALSE,
    showticklabels = TRUE,
    autotick = FALSE,
    ticks = "outside",
    tick0 = 0,
    dtick = 0.5,
    tickwidth = 1,
    position = -1,
    rangemode = "nonnegative"
  )
  
  p1 <-
    plot_ly(
      dat,
      x = ~ logFC,
      y = ~ NegLogPvalMod,
      type = "scatter",
      mode = "markers",
      color = ~ categ_Mod,
      size = ~ NegLogPvalMod,
      colors = c('#0C4B8E', '#BF382A'),
      hoverinfo = 'text',
      text = ~ paste(
        "Gene:",
        dat$gene,
        "</br>Fold Change:",
        logFC,
        "</br>-log 10[p-value]:",
        NegLogPvalMod
      )
    ) %>%
    layout(
      shapes = list(
        list(
          type = 'line',
          x0 = FC_Cutoff,
          x1 = FC_Cutoff,
          y0 = 0,
          y1 = ceiling(max(dat$NegLogPvalMod)),
          # y1 = 7.5, # change here for the heighest -log10(pval)
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
          x0 = -FC_Cutoff,
          x1 = -FC_Cutoff,
          y0 = 0,
          y1 = ceiling(max(dat$NegLogPvalMod)),
          # y1 = 7.5, # change here for the heighest -log10(pval)
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
          # x0 = -10, # change here for the heighest log2(fold change)
          # x1 = 10,  # change here for the heighest log2(fold change)
          x0 = -ceiling(max(abs(min(dat$logFC)), max(dat$logFC))),
          x1 = ceiling(max(abs(min(dat$logFC)), max(dat$logFC))),
          y0 = -log10(0.05),
          y1 = -log10(0.05),
          line = list(dash = 'dot', width = 1)
        )
      ),
      title = paste("volcano plot of moderated p-values"),
      xaxis = x_axis,
      yaxis = y_axis,
      showlegend = TRUE
    )
  
  p2 <-
    plot_ly(
      dat,
      x = ~ logFC,
      y = ~ NegLogPvalOrd,
      type = "scatter",
      mode = "markers",
      color = ~ categ_Ord,
      size = ~ NegLogPvalOrd,
      colors = c('#0C4B8E', '#BF382A'),
      hoverinfo = 'text',
      text = ~ paste(
        "Gene:",
        dat$gene,
        "</br>Fold Change:",
        logFC,
        "</br>-log 10[p-value]:",
        NegLogPvalOrd
      )
    ) %>%
    layout(
      shapes = list(
        list(
          type = 'line',
          x0 = FC_Cutoff,
          x1 = FC_Cutoff,
          y0 = 0,
          y1 = ceiling(max(dat$NegLogPvalOrd)),
          # y1 = 7.5, # change here for the heighest -log10(pval)
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
          x0 = -FC_Cutoff,
          x1 = -FC_Cutoff,
          y0 = 0,
          y1 = ceiling(max(dat$NegLogPvalOrd)),
          # y1 = 7.5, # change here for the heighest -log10(pval)
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
          # x0 = -10, # change here for the heighest log2(fold change)
          # x1 = 10, # change here for the heighest log2
          x0 = -ceiling(max(abs(min(dat$logFC)), max(dat$logFC))),
          x1 = ceiling(max(abs(min(dat$logFC)), max(dat$logFC))),
          y0 = -log10(0.05),
          y1 = -log10(0.05),
          line = list(dash = 'dot', width = 1)
        )
      ),
      title = paste("volcano plot of ordinary p-values"),
      xaxis = x_axis,
      yaxis = y_axis,
      showlegend = TRUE
    )
  
  filename_mod <- paste(filename_mod, '.html', sep = '')
  filename_ord <- paste(filename_ord, '.html', sep = '')
  htmlwidgets::saveWidget(as_widget(p1), filename_mod)
  htmlwidgets::saveWidget(as_widget(p2), filename_ord)
  
  # htmlwidgets::saveWidget(as_widget(p1), "Volcano_Limma_Mod_Pval.html")
  # htmlwidgets::saveWidget(as_widget(p2), "Volcano_Limma_Ord_Pval.html")
}
###########################################################################################################
###########################################################################################################
###########################################################################################################