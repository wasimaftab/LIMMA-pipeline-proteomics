median_normalization <- function(X, spike_in_rows = NULL){
  if(is.matrix(X)){
    stopifnot(length(dim(X)) == 2)
    if(is.null(spike_in_rows)){
      # Use all rows for normalization
      spike_in_rows <- seq_len(nrow(X))
    }else if(! is.numeric(spike_in_rows) && ! is.logical(spike_in_rows)){
      stop("spike_in_rows must either be a numeric or a logical vector")
    }else if(is.logical(spike_in_rows)){
      if(length(spike_in_rows) != nrow(X)){
        stop("The spike_in_rows logical vector must have one entry for each row.")
      }else if(all(spike_in_rows == FALSE)){
        stop("Not all elements of the spike_in_rows vector can be FALSE.")
      }
    }
    Xnorm <- X
    for(idx in seq_len(ncol(X))){
      Xnorm[, idx] <- X[, idx, drop=FALSE] -
        median(X[spike_in_rows, idx, drop=FALSE] - rowMeans(X[spike_in_rows, , drop=FALSE], na.rm=TRUE), na.rm=TRUE )
    }
    Xnorm
  }
}
################################################################################
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
################################################################################
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
################################################################################
# read_k_out_of_N <- function(rep_treats, group)
# {
#   str <- paste0("Enter an integer greater than 0 but less than the number of ", group," replicates(", rep_treats, ") = ")
#   n <- as.integer(readline(prompt = str))
#   if (is.na(n) || (n <= 0)) {
#     print("Enter positive integer only")
#     n <- readinteger(str)
#   } else if(n > rep_treats){
#     print(paste("Enter positive integer less than", rep_treats))
#     n <- readinteger(rep_treats)
#   }
#   return(n)
# }
###
## Sanity check k out of N replicates
get_k_out_of_N <- function(N, exp_type, type_molecule = "protein") {
  while (1) {
    k <-
      as.integer(readline(
        paste(
          'Enter k, to select only',
          paste0(type_molecule, "s"),
          'that have k non zero values out of',
          N,
          exp_type,
          'replicates and k>0 = '
        )
      ))
    if (k > N) {
      cat('k must be an integer>0 and <=', N, 'TRY AGAIN\n')
    } else if (k <= 0) {
      cat('k must be an integer>0 and <=', N, 'TRY AGAIN\n')
    } else {
      break
    }
  }
  return(k)
}
################################################################################
readnumber <- function(str)
{
  n <- readline(prompt = str)
  n <- as.integer(n)
  if (is.na(n)) {
    print("This variable expects an integer only")
    n <- readnumber(str)
  }
  return(n)
}
################################################################################
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
################################################################################
readfloat_0_1 <- function(str)
{
  n <- readline(prompt = str)
  n <- as.double(n)
  if (is.na(n) || (n <= 0)) {
    writeLines("Enter a number>0 & <=1, try again")
    n <- readfloat(str)
  }else if (n > 1){
    writeLines("Enter a number>0 & <=1, try again")
  }
  return(n)
}
################################################################################
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
################################################################################
data_sanity_check <- function(temp, exprement, exp_str) {
  # browser()
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
################################################################################
display_plotly_figs <-
  function(dat,
           FC_Cutoff,
           filename_mod,
           filename_ord) {
    m <- list(l = 60,
              r = 40,
              b = 55,
              t = 40)

    f <- list(family = "Arial, sans-serif",
              size = 20,
              # color = "#7f7f7f")
              color = 'black')

    f2 <- list(family = "Arial, sans-serif",
               size = 14,
               color = 'black')

    x_axis <- list(
      # constraintoward = "bottom",
      zeroline = TRUE,
      showline = FALSE,
      zerolinecolor = toRGB("black"),
      zerolinewidth = 2,
      title = "<b>log2 fold change</b>",
      # title = paste0(c(rep("&nbsp;", 2),"<b>log2 fold change</b>")),
      titlefont = f,
      showgrid = FALSE,
      showticklabels = TRUE,
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 1,
      tickwidth = 1,
      position = -1,
      rangemode = "tozero",
      tickfont = f2
    )

    y_axis <- list(
      zeroline = TRUE,
      showline = FALSE,
      zerolinecolor = toRGB("black"),
      zerolinewidth = 2,
      title = "<b>-log10 p-value</b>",
      titlefont = f,
      showgrid = FALSE,
      showticklabels = TRUE,
      autotick = FALSE,
      ticks = "outside",
      tick0 = 0,
      dtick = 0.5,
      tickwidth = 1,
      position = -1,
      rangemode = "nonnegative",
      tickfont = f2
    )

    # Plot Limma result
    p1 <-
      plot_ly(
        dat,
        x = ~ logFC,
        y = ~ NegLogPvalMod,
        type = "scatter",
        mode = "markers",
        color = ~ categ_Mod,
        colors = c('#0C4B8E', '#BF382A'),
        size = ~ abs(logFC),
        hoverinfo = 'text',
        showlegend = FALSE,
        text = ~ paste(
          "Uniprot:",
          dat$Uniprot,
          # Changed here from gene to Uniprot
          "</br>",
          "</br>Gene:",
          dat$Symbol,
          # Added here ProteinNames
          "</br>",
          "</br>Fold Change:",
          logFC,
          "</br>-log 10[p-value]:",
          NegLogPvalMod
        )
      ) %>%
      add_trace(marker = list(line = list(
        color = toRGB('black'),
        width = 1
      )),
      showlegend = TRUE) %>%
      layout(
        shapes = list(
          list(
            type = 'line',
            x0 = FC_Cutoff,
            x1 = FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalMod)),
            width = 2,
            # line = list(dash = 'dot', width = 1)
            line = list(
              dash = 'dot',
              width = 2,
              color = 'black'
            )
          ),
          list(
            type = 'line',
            x0 = -FC_Cutoff,
            x1 = -FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalMod)),
            # line = list(dash = 'dot', width = 1)
            line = list(
              dash = 'dot',
              width = 2,
              color = 'black'
            )
          ),
          list(
            type = 'line',
            x0 = -ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            x1 = ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            y0 = -log10(0.05),
            y1 = -log10(0.05),
            # line = list(dash = 'dot', width = 1)
            line = list(
              dash = 'dot',
              width = 2,
              color = 'black'
            )
          )
        ),
        title = paste("volcano plot of moderated p-values"),
        xaxis = x_axis,
        yaxis = y_axis,
        showlegend = TRUE,
        titlefont = f,
        margin = m
      )


    # Plot Ordinary t-test result
    p2 <-
      plot_ly(
        dat,
        x = ~ logFC,
        y = ~ NegLogPvalOrd,
        type = "scatter",
        mode = "markers",
        color = ~ categ_Ord,
        colors = c('#0C4B8E', '#BF382A'),
        size = ~ abs(logFC),
        hoverinfo = 'text',
        showlegend = FALSE,
        text = ~ paste(
          "Uniprot:",
          dat$Uniprot,
          # Changed here from gene to Uniprot
          "</br>",
          "</br>Gene:",
          dat$Symbol,
          # Added here ProteinNames
          "</br>",
          "</br>Fold Change:",
          logFC,
          "</br>-log 10[p-value]:",
          NegLogPvalMod
        )
      ) %>%
      add_trace(marker = list(line = list(
        color = toRGB('black'),
        width = 1
      )),
      showlegend = TRUE) %>%
      layout(
        shapes = list(
          list(
            type = 'line',
            x0 = FC_Cutoff,
            x1 = FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalMod)),
            width = 2,
            # line = list(dash = 'dot', width = 1)
            line = list(
              dash = 'dot',
              width = 2,
              color = 'black'
            )
          ),
          list(
            type = 'line',
            x0 = -FC_Cutoff,
            x1 = -FC_Cutoff,
            y0 = 0,
            y1 = ceiling(max(dat$NegLogPvalMod)),
            # line = list(dash = 'dot', width = 1)
            line = list(
              dash = 'dot',
              width = 2,
              color = 'black'
            )
          ),
          list(
            type = 'line',
            x0 = -ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            x1 = ceiling(max(abs(min(
              dat$logFC
            )), max(dat$logFC))),
            y0 = -log10(0.05),
            y1 = -log10(0.05),
            # line = list(dash = 'dot', width = 1)
            line = list(
              dash = 'dot',
              width = 2,
              color = 'black'
            )
          )
        ),
        title = paste("volcano plot of ordinary p-values"),
        xaxis = x_axis,
        yaxis = y_axis,
        # showlegend = TRUE,
        titlefont = f,
        margin = m
      )

    # subDir <- "Results"
    subDir <- paste0("Results_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    mainDir <- getwd()

    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))

    # ## Save as .html
    htmlwidgets::saveWidget(as_widget(p1), paste(filename_mod, '.html', sep = ''))
    htmlwidgets::saveWidget(as_widget(p2), paste(filename_ord, '.html', sep = ''))

    # ## Save as .pdf
    # plotly::orca(p1, file = paste(filename_mod, '.pdf', sep = ''))
    # plotly::orca(p2, file = paste(filename_ord, '.pdf', sep = ''))
  }
###########################################################################################################
