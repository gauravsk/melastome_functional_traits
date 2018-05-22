source("scripts/functions/add_trendline.R")

tiff("fig1.tiff",width=11,height=7,res=600,units="in",compression="lzw",bg="white",type="cairo")


par(mfrow = c(2,4), oma = c(8,1,1,1), mar = c(0, 5, 2, 1), 
    cex.lab = 1.5, cex.main = 2, xpd = F, cex.axis = 1.25)

lcex = 1.3
sp_by_location$pchvec = pchvec
sp_by_location$colvec = NA
cols = grey(alpha=0.7,c(1,.5,0))
for(i in 1:3){sp_by_location[sp_by_location$pchvec==(21:23)[i],"colvec"] = cols[i]}
head(sp_by_location)

# SLA
plot(sp_by_location$log_sla~jitter(sp_by_location$elevation, factor = .5),
     xlab = "", pch = sp_by_location$pchvec, bty = "l", xaxt = "n", 
     bg = sp_by_location$colvec, cex = 1.5, ylim = c(1.9,2.9),
     ylab = expression(paste('log'[10],'(SLA) (','cm'^2,'g'^-1,')')), las = 1)
add_line(trait = "log_sla", df = sp_by_location, log = F)
#box(lwd = 0.5)
mtext(side = 3, text = "a", adj = 0.02, line = -1, font = 2, cex = lcex)

# LDMC
plot(sp_by_location$log_ldmc~jitter(sp_by_location$elevation, factor = .5), 
     xlab = "", pch = sp_by_location$pchvec, bty = "l", xaxt = "n", 
     bg = sp_by_location$colvec, cex = 1.5, ylim = c(1.05,1.75),
     ylab = expression(paste('log'[10],'(LDMC) (','mg g'^-1,')')), las = 1)
add_line(trait = "log_ldmc", df = sp_by_location, log = F)
#box(lwd = 0.5)
mtext(side = 3, text = "b", adj = 0.02, line = -1, font = 2, cex = lcex)

# Leaf Area
plot(sp_by_location$log_leaf_area~jitter(sp_by_location$elevation, factor = .5), 
     xlab = "", pch = sp_by_location$pchvec, bty = "l", xaxt = "n", 
     bg = sp_by_location$colvec, cex = 1.5,
     ylab = expression(paste('log'[10],'(Leaf area) (','cm'^2,')')), las = 1,
     ylim=c(.2,3.2))
add_line(trait = "log_leaf_area", df = sp_by_location, log = F)
#box(lwd = 0.5)
mtext(side = 3, text = "c", adj = 0.02, line = -1, font = 2, cex = lcex)

# Leaf force to punch
plot(sp_by_location$log_specific_force~jitter(sp_by_location$elevation, factor = .5), cex = 1.5,
     xlab = "", pch = sp_by_location$pchvec, bty = "l", 
     bg = sp_by_location$colvec,  las = 1, ylim = c(-0.2,0.9),
     ylab = expression(paste('log'[10],'(Leaf force to punch) (','N m'^-2,')')))
add_line(trait = "log_specific_force", df = sp_by_location, log = F, pos=0.1)
#box(lwd = 0.5)
mtext(side = 3, text = "d", adj = 0.02, line = -1, font = 2, cex = lcex)

# Leaf N
plot(sp_by_location$lnc~jitter(sp_by_location$elevation, factor = .5), 
     xlab = "", pch = sp_by_location$pchvec, bty = "l", 
     bg = sp_by_location$colvec, cex = 1.5, las = 1, ylab = "Leaf nitrogen concentration (%)")
add_line(trait = "lnc", df = sp_by_location, log = F)
#box(lwd = 0.5)
mtext(side = 3, text = "e", adj = 0.02, line = -1, font = 2, cex = lcex)

# Stem Density
plot(sp_by_location$stem_density~jitter(sp_by_location$elevation, factor = .5),
     xlab = "", pch = sp_by_location$pchvec, bty = "l", 
     bg = sp_by_location$colvec,  las = 1, cex = 1.5,
     ylab = expression(paste('Stem-specific density (','g cm'^-3,')')))
add_line(trait = "stem_density", df = sp_by_location, log = F, pos=.65)
#box(lwd = 0.5)
mtext(side = 3, text = "f", adj = 0.02, line = -1, font = 2, cex = lcex)

# Seed size
plot(sp_by_location$log_seed_area~jitter(sp_by_location$elevation, factor = .5), 
     xlab = "", pch = sp_by_location$pchvec, bty = "l", bg = sp_by_location$colvec, cex = 1.5, las = 1,
     ylab = expression(paste('log'[10],'(Seed area) (',mu,'m'^2,')')), 
     ylim = c(4.5,8))
add_line(trait = "log_seed_area", df = sp_by_location, log = F)
#box(lwd = 0.5)
mtext(side = 3, text = "g", adj = 0.02, line = -1, font = 2, cex = lcex)

# Legend
plot(sp_by_location$log_seed_area~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", bty = "n", type="n", ylab = "", xaxt="n", yaxt="n")
legend(x=0, y=6.5, pch = c(22,23,21), 
       pt.bg = cols[c(2,3,1)], xpd = NA, y.intersp = 1,
       legend = levels(sp_by_location$growth_form), horiz = F, 
       bty = "n", pt.cex = 1.5, cex = 1.7, x.intersp = 0.8)

# Title and axis
# title("Trait turnover across elevation", outer = T, line = 1, cex = 1.25)
mtext(side = 1, text = "Altitude (m asl)", outer = T, line = 5, cex = 1.25)


dev.off()
