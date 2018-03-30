# Making Figure 1 for melastome draft
## Make a conceptual figure...
library(plotrix)
make_pdf <- FALSE
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}
# cols <- c("steelblue", "#FFBB00")
# cols <- c("#75B1A9", "#D9B44A")
cols <- c("#4897D8", "#F8A055")
cols <- add.alpha(cols, alpha = 0.75)

cols <- c("white", "gray60")

if(make_pdf == TRUE) {pdf("~/Desktop/test1.pdf", width = 8, height = 5)}
# par(mfrow = c(2,2), mar = c(1, 1, 1, 1), oma = c(2.3, 2.3, 3, 0))
par(mfrow = c(1,3), mar = c(1, 1, 1, 1), oma = c(2.3, 2.3, 3, 0))

# Panel A
# plot(1, type = "n", xaxt = "n", yaxt = "n")
# draw.ellipse(x = 1, y = 1, .18,.41, angle = 90, col = cols[1])
# draw.ellipse(x = 1, y = 1, .08,.4, angle = 90, col = cols[2])
# text(.61, 1.36, "A", cex = 1.5)
# legend("topright", bty = "n", pch = 21, col = "black", pt.bg = cols, pt.cex = 1.25, legend = c("Community", "Focal Clade"))

# Panel B
plot(1, type = "n", xaxt = "n", yaxt = "n")
draw.ellipse(x = 1, y = 1, .22,.5, angle = -45, col = cols[1])
draw.ellipse(x = 1, y = 1, .1,.24, angle = 90, col = cols[2])
# text(.61, 1.36, "B", cex = 1.5)

# Panel C
plot(1, type = "n", xaxt = "n", yaxt = "n")
draw.ellipse(x = 1, y = 1, .22,.5, angle = -45, col = cols[1])
draw.ellipse(x = 1, y = 1, .08,.45, angle = -45, col = cols[2])
# text(.61, 1.36, "C", cex = 1.5)

# Panel D
plot(1, type = "n", xaxt = "n", yaxt = "n")
draw.ellipse(x = 1, y = 1, .22,.5, angle = -45, col = cols[1])
draw.ellipse(x = 1, y = 1, .08,.2, angle = 45, col = cols[2])
# text(.61, 1.36, "D", cex = 1.5)

mtext(outer = T, side = 1, text = "Elevation", cex = 1.25)
mtext(outer = T, side = 2, text = "Trait value", cex = 1.25)
title("Schematic of community- and clade-wide trait turnover across gradients", outer = "T", cex = 1.5)
if(make_pdf == TRUE){dev.off()}