# phy_plottr
# used to make figure of phylogeny with traits
# Gaurav Kandlikar, gkandlikar@ucla.edu

phy_plottr <- function (df, trait, phy, title, log = T, checkorder = T, pch = 19, cex = 1,...) {
  
  
  if (checkorder) {
    df <- df[(phy$tip.label[phy$edge[phy$edge[,2]<= Ntip(phy),2]]),]
  }
  ntips <- length(phy$tip.label)
  traitcol <- df[, which(colnames(df) == trait)]
  if (log) {traitcol <- log10(traitcol)}
  names(traitcol) <- rownames(df)
  
  
  to_plot <- data.frame(x = traitcol, y = 2:(ntips+1))
  plot(x = na.omit(to_plot$x), y = to_plot[which(!is.nan(to_plot$x)), "y"], pch = pch, cex = cex, axes = F, type = "b", xlab = "", ylab = "", ylim = c(1,ntips))
  axis(side = 1)
  mtext(side = 3, line = 1, text = title, cex = .8, ...)
  
  # Update the numbers in this after talking to Marcel since he came up with the numbers
  # for (i in c(1, 11, 20, 24, 28, 43, 61, 64, 70, 85, 92)) {
  #   abline(h = i+.5, lty = 1, lwd = 1, col = "grey50")
  # }
  
  # Plot the grid
  # for (i in 1:length(phy$tip.label)) {
  #   if (i %% 5 == 0) {abline (h = i, lwd = 0.75, col = "grey50")}
  #   else {abline(h = i, lty = 3, lwd = 0.75, col = "grey50")}
  # }
  

  #################################################
  # This was the old way I did things:----------- #
  #################################################
  
  # plot(1, type = "n", ylim = c(0, ntips-1),
  #      xlim = c(min(traitcol, na.rm = T) - .25,max(traitcol, na.rm = T) + .25),
  #      axes = F, ylab = "", xlab = "")
  # axis(side = 1)
  # for (i in 2:length(phy$tip.label)) {
  #   points(x = na.omit(c(traitcol[i], traitcol[i-1])),
  #          y = c(i, i-1), type = "o", pch = pch, cex = cex, ...)
  # }

}