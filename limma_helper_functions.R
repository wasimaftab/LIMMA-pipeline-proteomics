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

readinteger_binary <- function(str)
{
  n <- readline(prompt = str)
  n <- as.integer(n)
  if (is.na(n) || (n < 0) || (n > 1)) {
    print("Enter positive integer only")
    n <- readinteger_binary(str)
  }
  return(n)
}

data_sanity_check <- function(temp, exprement, exp_str) {
  exprement_reps <-
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
      colors = c('#0C4B8E', '#BF382A'),
      hoverinfo = 'text',
      text = ~ paste(
        "Gene:",
        dat$gene,
        "</br>",
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
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
          x0 = -FC_Cutoff,
          x1 = -FC_Cutoff,
          y0 = 0,
          y1 = ceiling(max(dat$NegLogPvalMod)),
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
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
      colors = c('#0C4B8E', '#BF382A'),
      hoverinfo = 'text',
      text = ~ paste(
        "Gene:",
        dat$gene,
        "</br>",
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
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
          x0 = -FC_Cutoff,
          x1 = -FC_Cutoff,
          y0 = 0,
          y1 = ceiling(max(dat$NegLogPvalOrd)),
          line = list(dash = 'dot', width = 1)
        ),
        list(
          type = 'line',
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

  subDir <- "Results"
  mainDir <- getwd()
  
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  # ## Save as .html
  htmlwidgets::saveWidget(as_widget(p1), paste(filename_mod, '.html', sep = ''))
  htmlwidgets::saveWidget(as_widget(p2), paste(filename_ord, '.html', sep = ''))
  
  # ## Save as .pdf
  export(p1, file = paste(filename_mod, '.pdf', sep = ''))
  export(p2, file = paste(filename_ord, '.pdf', sep = ''))

}
###########################################################################################################
###########################################################################################################
###########################################################################################################