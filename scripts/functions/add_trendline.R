# A function to add a trendline to a plot if a linear regression is significant
# df is the data frame containing trait data
#     NOTE: this function also needs an "elevation" column in df
# trait is the name of the column containing trait data
# log = T/F based on whether the trait should be log10 transformed prior to lm()
# digits = number of digits to show in R and p

add_line <- function(trait, df, log = F, digits = 3) {
  if (log == F) {current_trt <- unlist(df[,trait])}
  else {current_trt <- log10(unlist(df[,trait]))}
  net_lm <- lm(current_trt~df$elevation)
  rsq <- round(summary(net_lm)$r.squared, 2)
  pval <- anova(net_lm)$'Pr(>F)'[1]
  
  if (pval < 0.05) {
    abline(net_lm, lwd = 1.5)

    to_write <- bquote(atop(italic(R)^2 == .(rsq), 
                        italic(p) == .(format(pval, digits = 2))))
    text(x = 1500, y = par("usr")[4], pos = 1, offset = 1, labels = to_write, bty = 'n', cex = 1.5)
  }
}

