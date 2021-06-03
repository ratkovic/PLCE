#' Assessing Positivity of a Treatment Variable
#'
#' @param obj A PLCE object.
#' @param trim The percentage to trim off each end of the excess kurtosis statistic.
#' This helps in reducing
#' @export
#'
#'
#' @return \describe{
#' \item{diff}{The estimated excess kurtosis per datum.}
#'}
#'
#' @references Ratkovic, Marc.  2021. "Relaxing Assumptions, Improving Inference:
#' Utilizing Machine Learning for Valid Causal Inference."  Working Paper.
#'
#' @rdname pos_plot


pos_plot <- function(obj, trim = 0.025) {
  p1 <- pos_measure(obj, trim)
  
  
  
  gg_color <- function(n, alpha0 = 1) {
    hues = seq(15, 375, length = n + 1)
    hcl(
      h = hues,
      l = 65,
      c = 100,
      alpha = alpha0
    )[1:n]
  }
  
  
  
  cols <- gg_color(2)
  
  par(
    mar = c(3, 3, 1, 0.2),
    # Dist' from plot to side of page
    mgp = c(2, 0.4, 0),
    # Dist' plot to label
    las = 1,
    # Rotate y-axis text
    tck = -.01,
    # Reduce tick length
    xaxs = "i",
    yaxs = "i",
    # Remove plot padding
    oma = c(.2, .2, 0, 0)
  )
  
  
  r1 <- diff(range(c(p1$diff, 0)))
  plot(
    make.rank(p1$diff),
    p1$diff,
    type = "n",
    col = "red",
    lwd = 2,
    ylim = range(c(p1$diff, 0)),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    xlim = c(-5, 105)
  )
  #abline(c(0,1),col="gray50")
  abline(h = 0, lwd = 2)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "gray90")
  abline(h = (-100:100) / 10, col = "white")
  abline(v = (0:100) * 10, col = "white")
  
  mtext(
    side = 1,
    line = 1.7,
    font = 2,
    text = "Percentile",
    las = 0
  )
  mtext(
    side = 2,
    line = 2,
    font = 2,
    text = "Excess Kurtosis",
    las = 0
  )
  axis(1)
  axis(2)
  lines(make.rank(p1$diff),
        p1$diff,
        col = cols[1],
        lwd = 2)
  
  output <- list("ExcessKurtosis" = p1$diff)
  invisible(output)
}
