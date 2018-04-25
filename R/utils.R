scoreIPSCustom <- function(x, IPS.custom, IPS20 = TRUE){
  d <- x[x > 0]
  names(d) <- names(x)[x > 0]
  ips <- IPS.custom[names(d),]
  tab <- data.frame(ips, d)
  avs <- apply(tab, 1, function(x) prod(x))
  av <- apply(tab, 1, function(x) x[2] * x[3])
  res <- sum(avs, na.rm = T) / sum(av, na.rm = T)
  if(IPS20){
    res <- (4.75 * res) - 3.75
  }
  return(res)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  abline(a = 0, b = 1, lty = "dashed")
  points(x, y, pch = 20)
}

panel.txt <- function(x, y, digits = 2, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  r <- format(c(r, 0.123456789), digits = digits)[1]
  rpval <- cor.test(x, y)$p.value
  rpval <- ifelse(rpval < 10^-5, "*", "")
  mse <- sum((x - y)^2)/length(x)
  mse <- format(c(mse, 0.123456789), digits = digits)[1]
  txt <- paste0("cor = ", r, rpval, "\n", "MSE= ", mse)
  text(0.5, 0.5, txt)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 2) )
  h <- hist(x, plot = FALSE, breaks = seq(0, 5, 0.25))
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$density; #y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey60", ...)
}