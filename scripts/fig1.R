source("scripts/functions/add_trendline.R")

tiff("fig1.tiff",width=11,height=7,res=600,units="in",compression="lzw",bg="white",type="cairo")


par(mfrow = c(2,4), oma = c(8,1,5,1), mar = c(0, 5, 0, 1), 
    cex.lab = 1.5, cex.main = 2, xpd = F, cex.axis = 1.25)

# SLA
plot(sp_by_location$log_sla~jitter(sp_by_location$elevation, factor = .3),
     xlab = "", pch = pchvec, bty = "n", xaxt = "n", bg = colvec, cex = 1.5,
     ylab = expression(paste('log'[10],'(SLA) (','cm'^2,'g'^-1,')')),
     ylim=c(1.9,2.9))
add_line(trait = "log_sla", df = sp_by_location, log = F)
box(lwd = 0.5)
#mtext(side = 3, text = "A", adj = 0.02, line = -1.3, font = 2)

# LDMC
plot(sp_by_location$log_ldmc~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", pch = pchvec, bty = "n", xaxt = "n", bg = colvec, cex = 1.5,
     ylab = expression(paste('log'[10],'(LDMC) (','mg g'^-1,')')),
     ylim=c(1.05,1.75))
add_line(trait = "log_ldmc", df = sp_by_location, log = F)
box(lwd = 0.5)
#mtext(side = 3, text = "B", adj = 0.02, line = -1.3, font = 2)

# Leaf Area
plot(sp_by_location$log_leaf_area~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", pch = pchvec, bty = "n", xaxt = "n", bg = colvec, cex = 1.5,
     ylab = expression(paste('log'[10],'(Leaf Area) (','cm'^2,')')),
     ylim=c(.2,3.1))
add_line(trait = "log_leaf_area", df = sp_by_location, log = F)
box(lwd = 0.5)
#mtext(side = 3, text = "C", adj = 0.02, line = -1.3, font = 2)

# Specific Force
plot(sp_by_location$log_specific_force~jitter(sp_by_location$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = pchvec, bty = "n", bg = colvec, 
     ylab = expression(paste('log'[10],'(Leaf Toughness) (','N m'^-2,')')))
add_line(trait = "log_specific_force", df = sp_by_location, log = F, pos=0.1)
box(lwd = 0.5)
#mtext(side = 3, text = "D", adj = 0.02, line = -1.3, font = 2)

# Leaf N
plot(sp_by_location$lnc~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", pch = pchvec, bty = "n", bg = colvec, cex = 1.5,
     ylab = "Leaf Nitrogen (%)")
add_line(trait = "lnc", df = sp_by_location, log = F)
box(lwd = 0.5)
#mtext(side = 3, text = "E", adj = 0.02, line = -1.3, font = 2)

# Stem Density
plot(sp_by_location$stem_density~jitter(sp_by_location$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = pchvec, bty = "n", bg = colvec, 
     ylab = expression(paste('Stem density (','g cm'^-3,')')))
add_line(trait = "stem_density", df = sp_by_location, log = F, pos=.65)
box(lwd = 0.5)

# Seed size
plot(sp_by_location$log_seed_area~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", pch = pchvec, bty = "n", bg = colvec, cex = 1.5,
     ylab = expression(paste('log'[10],'(Seed area) (',mu,'m'^2,')')))
add_line(trait = "log_seed_area", df = sp_by_location, log = F)
box(lwd = 0.5)

# Legend
plot(sp_by_location$log_seed_area~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", bty = "n", bg = colvec, type="n", ylab = "", xaxt="n", yaxt="n")

legend(x=0, y=6.5, pch = pchs, 
       pt.bg = colors_notransp, xpd = NA, y.intersp = 1,
       legend = levels(sp_by_location$growth_form), horiz = F, 
       bty = "n", pt.cex = 1.5, cex = 1.7, x.intersp = 0.8)
#mtext(side = 3, text = "F", adj = 0.02, line = -1.3, font = 2)

# Title and axis
# title("Trait turnover across elevation", outer = T, line = 1, cex = 1.25)
mtext(side = 1, text = "Elevation (m)", outer = T, line = 5, cex = 1.25)


dev.off()
