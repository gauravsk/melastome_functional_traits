# Complete analysis of the melastome trait turnover data
# Gaurav Kandlikar, gkandlikar@ucla.edu
# Last updated 12 April 2018

# Set up ----------
library(dplyr); library(ape); library(picante); library(vegan)
# setwd("~/grad/Dropbox/melastome_trait_data/biotropica_submission/")
individual_traits <- read.csv("data/individual_level_traits.csv")
tbl_df(individual_traits)
default_pars <- par()
# Toggle these to change some downstream things
woody_only <- FALSE
reps_for_K <- 999
write_outputs <- FALSE

# Data munging -------
# Summarize the traits into a species-by-elevation table: 
sp_by_location <- 
  individual_traits %>% 
  group_by(full_name, genus, species, location, elevation, growth_form) %>% 
  summarise(sla = mean(sla, na.rm = T), 
            ldmc = mean(ldmc, na.rm = T), 
            leaf_area = mean(area_actual, na.rm = T), 
            specific_force = mean(fracture, na.rm = T), 
            stem_density = mean(stem_density, na.rm = T), 
            thickness = mean(thickness, na.rm = T))

# Read in and process Leaf N data
lncs <- read.csv("data/leaf_N.csv"); nrow(lncs)

# Make a new column for a future species X site merge
sp_by_location$full_name_site <- paste(sp_by_location$full_name, sp_by_location$location)
lncs$full_name_site <- paste(lncs$species, lncs$Location)

# Apparently there are two more rows in the lnc dataframe than in sp_by_location
# Figure out which those are 
N_but_no_traits <- which(is.na(match(lncs$full_name_site, sp_by_location$full_name_site)))
# These are the species X sites for which we lack trait data
lncs[N_but_no_traits, c("species", "Location")]
# we lacked trait data for those two- so let's drop them.
lncs <- lncs[-N_but_no_traits,]
tbl_df(lncs)


# Subset the dataset to only include Mean N Percent and the row on which to merge
# Also, rename "Mean.N_Percent" to just "lnc"
lncs <- lncs %>% select(lnc = Mean.N_Percent., full_name_site)
sp_by_location <- merge(sp_by_location, lncs, by = "full_name_site")
# Drop the full_name_site column, which contain redundant info
sp_by_location <- sp_by_location %>% select(-full_name_site)
tbl_df(sp_by_location)
str(sp_by_location)


# Include area-based leaf nitrogen content (mg/cm2)
sp_by_location$alnc <- sp_by_location$lnc/(100*sp_by_location$sla)


# Read in Seed size data 
seed <- read.csv("data/seed_size.csv")
tbl_df(seed)
plot(seed$mean_width~seed$mean_length)
# Since mean length and mean width are two different columns right now, 
# I will combine them into a "seed area"
# To do this I model the seed as a rectangle...
seed$seed_area <- seed$mean_width*seed$mean_length
# Let's confirm that this correlates well with the measured dimensions
par(mfrow = c(1,2))
plot(seed$seed_area~seed$mean_width)
plot(seed$seed_area~seed$mean_length)
# The two are non-linearly correlated- which makes sense, since we're multiplying...
# I will do the analysis of seed size as a projected rectangle, but should confirm with NK and RK and FM...


# Lets subset the dataframe to just the seed area- i.e. drop out redundant width/length
# We are also losing info like fruit type, which is not used in the analysis anywhere.
seed <- seed %>% select(full_name, seed_area)
# Merge seed dataframe with the rest of the traits
sp_by_location <- left_join(sp_by_location, seed, by = "full_name")
# NOTE!
# Since Seed data is only at the per species level (and not at the Species X Site level), 
# species values for this trait will be equal at all elevations at which the species occurs

# For example....
sp_by_location %>% 
  group_by(full_name) %>% filter(n()>1) %>% select(full_name, elevation, sla, seed_area)


# OK! Take a look at the current sp_by_location dataset, since this will be the basis of analysis...
View(sp_by_location)

# If trees_only is TRUE, then subset the dataset by growth form
if (woody_only == TRUE) {
  sp_by_location <- sp_by_location %>% filter(growth_form == "woody")
}

# Read in the phylogeny
phylo <- read.tree("data/melas_cr_july2016.tre")
# Get rid of the outgroup
# Tips that are in the phylogeny but not in the dataset:
mismatch <- phylo$tip.label[which(!phylo$tip.label %in% sp_by_location$full_name)]
phylo <- drop.tip(phylo, mismatch)

# Begin trait analysis ----------
# First, check distribution of traits
par(mfrow = c(2,5))
mapply(hist, sp_by_location %>% select(sla, ldmc, leaf_area, 
                                       specific_force, stem_density, 
                                       thickness, lnc, alnc, seed_area), 
       main = colnames(sp_by_location %>% select(sla, ldmc, leaf_area, 
                                                 specific_force, stem_density, 
                                                 thickness, lnc, alnc, seed_area)))
# It looks like SLA, LDMC, leaf area, thickness, and specific force need to be log-transformed for normalization
sp_by_location <- sp_by_location %>% mutate(log_sla = log10(sla), 
                                            log_ldmc = log10(ldmc), 
                                            log_leaf_area = log10(leaf_area), 
                                            log_specific_force = log10(specific_force), 
                                            log_seed_area = log10(seed_area), 
                                            log_thickness = log10(thickness),
                                            log_alnc = log10(alnc))

par(mfrow = c(2,5))
mapply(hist, sp_by_location %>% select(log_sla, log_ldmc, log_leaf_area, 
                                       log_specific_force, stem_density, 
                                       log_thickness, lnc, log_alnc, log_seed_area), 
       main = colnames(sp_by_location %>% select(log_sla, log_ldmc, log_leaf_area, 
                                                 log_specific_force, stem_density, 
                                                 log_thickness, lnc, log_alnc, log_seed_area)))

# Make a trait turnover plot ----------
# First, source in the add_trendline() function 
source("scripts/functions/add_trendline.R")


# Plot setup
color <- grey(.5, alpha = 0.5)
pchs <- c(22, 23, 21)
pchvec <- pchs[sp_by_location$growth_form]
# colors_notransp <- c("#F96876", "#5F7CC6", "#FFF76B")
colors_notransp <- c("#696969")
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

colors <- add.alpha(colors_notransp, alpha = 0.5)
colvec <- colors
# colvec <- colors[(sp_by_location$growth_form)]

if (write_outputs == TRUE) {pdf("figures/trait_turnover.pdf", height = 15, width = 15)}

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
     ylim=c(1.05,1.7))
add_line(trait = "log_ldmc", df = sp_by_location, log = F)
box(lwd = 0.5)
#mtext(side = 3, text = "B", adj = 0.02, line = -1.3, font = 2)

# Leaf Area
plot(sp_by_location$log_leaf_area~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", pch = pchvec, bty = "n", xaxt = "n", bg = colvec, cex = 1.5,
     ylab = expression(paste('log'[10],'(Leaf Area) (','cm'^2,')')))
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

# Leaf N (area-based)
#plot(sp_by_location$log_alnc~jitter(sp_by_location$elevation, factor = .3), 
#     xlab = "", pch = pchvec, bty = "n", bg = colvec, cex = 1.5,
#     ylab = "Leaf Nitrogen (g/cm2)")
#add_line(trait = "log_alnc", df = sp_by_location, log = F)
#box(lwd = 0.5)
#mtext(side = 3, text = "E", adj = 0.02, line = -1.3, font = 2)

# Stem Density
plot(sp_by_location$stem_density~jitter(sp_by_location$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = pchvec, bty = "n", bg = colvec, 
     ylab = expression(paste('Stem density (','g cm'^-3,')')))
add_line(trait = "stem_density", df = sp_by_location, log = F, pos=.65)
box(lwd = 0.5)

# Seed size
tmp = subset(sp_by_location,log_seed_area<7)
plot(tmp$log_seed_area~jitter(tmp$elevation, factor = .3), 
     xlab = "", pch = pchvec, bty = "n", bg = colvec, cex = 1.5,
     ylab = expression(paste('log'[10],'(Seed area) (',mu,'m'^2,')')))
add_line(trait = "log_seed_area", df = tmp, log = F)
box(lwd = 0.5)

# Legend
plot(sp_by_location$log_seed_area~jitter(sp_by_location$elevation, factor = .3), 
     xlab = "", bty = "n", bg = colvec, type="n", ylab = "", xaxt="n", yaxt="n")

legend(x=0, y=6.5, pch = pchs, 
       pt.bg = colors_notransp, xpd = NA, y.intersp = 1,
       legend = levels(sp_by_location$growth_form), horiz = F, 
       bty = "n", pt.cex = 1.33, cex = 1.33, x.intersp = 0.8)
#mtext(side = 3, text = "F", adj = 0.02, line = -1.3, font = 2)

# Title and axis
# title("Trait turnover across elevation", outer = T, line = 1, cex = 1.25)
mtext(side = 1, text = "Elevation (m)", outer = T, line = 5, cex = 1.25)


if(write_outputs == TRUE) {dev.off()}
par(default_pars)


# # Projected Seed Area - not including this
 plot(sp_by_location$log_seed_area~jitter(sp_by_location$elevation, factor = .3), 
      xlab = "", pch = pchvec, bty = "n", xaxt = "n", bg = colvec, 
      ylab = expression(paste('log'[10],'(Seed area) (',mu,'m'^2,')')))
 add_line(trait = "log_seed_area", df = sp_by_location, log = F)
 box(lwd = 0.5)
 summary(lm(log_seed_area~elevation,sp_by_location))

 
# Trait turnover for woody species only --------
woody_only <- sp_by_location %>% filter(growth_form == "woody")
woody_only <- data.frame(woody_only)
if (write_outputs == TRUE) {pdf("figures/trait_turnover_woodyOnly.pdf", height = 15, width = 15)}
# Plot setup
par(mfrow = c(2,3), oma = c(5,1,5,0), mar = c(0, 5, 0, 1), 
    cex.lab = 1.5, cex.main = 2)
color <- grey(.5, alpha = 0.5)

# SLA
plot(woody_only$log_sla~jitter(woody_only$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = 21, bty = "n", xaxt = "n", bg = color, 
     ylab = expression(paste('log'[10],'(SLA) (','cm'^2,'g'^-1,')')))
add_line(trait = "log_sla", df = woody_only, log = F)
box(lwd = 0.5)

# LDMC
plot(woody_only$log_ldmc~jitter(woody_only$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = 21, bty = "n", xaxt = "n", bg = color, 
     ylab = expression(paste('log'[10],'(LDMC) (','mg g'^-1,')')))
add_line(trait = "log_ldmc", df = woody_only, log = F)
box(lwd = 0.5)

# Leaf Area
plot(woody_only$log_leaf_area~jitter(woody_only$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = 21, bty = "n", xaxt = "n", bg = color, 
     ylab = expression(paste('log'[10],'(Leaf Area) (','cm'^2,')')))
add_line(trait = "log_leaf_area", df = woody_only, log = F)
box(lwd = 0.5)

# Specific Force
plot(woody_only$log_specific_force~jitter(woody_only$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = 21, bty = "n", bg = color, 
     ylab = expression(paste('log'[10],'(Leaf Toughness) (','N m'^-2,')')))
add_line(trait = "log_specific_force", df = woody_only, log = F)
box(lwd = 0.5)

# Leaf N
plot(woody_only$lnc~jitter(woody_only$elevation, factor = .3), 
     xlab = "", pch = 21, bty = "n", bg = color, cex = 1.5,
     ylab = "Leaf Nitrogen (%)")
add_line(trait = "lnc", df = woody_only, log = F)
box(lwd = 0.5)

# Stem Density
plot(woody_only$stem_density~jitter(woody_only$elevation, factor = .3), cex = 1.5,
     xlab = "", pch = 21, bty = "n", bg = color, 
     ylab = expression(paste('Stem density (','g cm'^-3,')')))
add_line(trait = "stem_density", df = woody_only, log = F)
box(lwd = 0.5)

# Title and axis
title("Trait turnover across elevation (Woody species only)", outer = T, line = 1, cex = 1.25)
mtext(side = 1, text = "Elevation (m)", outer = T, line = 3, cex = 1.25)

if(write_outputs == TRUE) {dev.off()}
par(default_pars)

# # Projected Seed Area - not including this
plot(woody_only$log_seed_area~jitter(woody_only$elevation, factor = .3), 
     xlab = "", pch = pchvec, bty = "n", xaxt = "n", bg = colvec, 
     ylab = expression(paste('log'[10],'(Seed area) (',mu,'m'^2,')')))
add_line(trait = "log_seed_area", df = woody_only, log = F)
box(lwd = 0.5)
summary(lm(log_seed_area~elevation,woody_only))


# Mantel correlation of taxon/functional/phylo dissimilarities -------

# Elevation and Temperature dissimilarities

elev_diff <- dist(matrix(c(30, 500, 800, 2000, 2500)))
elev_diff
temp_diff <- dist(matrix(c(24.61551,22.20317,20.66338,14.50421,11.93789)))
temp_diff


# Taxonomic dissimilarity 
# First, make a species by site matrix
site_by_species <- reshape2::acast(sp_by_location, elevation~full_name)
site_by_species[,1:4]

# Make a list that shows which species are in each community
species_in_community <- list()
for (i in 1:5) {
  species_in_community[[i]] <- paste(names(site_by_species[i, site_by_species[i, ] > 0]), 
                                     rownames(site_by_species)[i], sep = "_")
}
# Number of species in each community
lapply(species_in_community, length)

# Compute the Binary Jaccard distance between each community

source("scripts/functions/gsk_jaccard.R")
tax_diff <- jaccard(site_by_species = site_by_species)

# Functional dissimilarity ----

# Summarize to species by elevation (i.e. combine "road" and "forest" at each site)
sp_by_elev <- sp_by_location %>% 
  group_by(full_name, genus, species, elevation, growth_form) %>% 
  summarise(log_sla = mean(log_sla, na.rm = T), 
            log_ldmc = mean(log_ldmc, na.rm = T), 
            log_leaf_area = mean(log_leaf_area, na.rm = T), 
            log_specific_force = mean(log_specific_force, na.rm = T), 
            stem_density = mean(stem_density, na.rm = T),
            lnc = mean(lnc), 
            log_seed_area = mean(log_seed_area, na.rm = T))

traits_by_elev <- sp_by_elev %>% ungroup () %>% select(log_sla:log_seed_area) 

rownames(traits_by_elev) <- paste(sp_by_elev$full_name, sp_by_elev$elevation, sep = "_")


# Make a separate vector for each trait to do separate trait turnover analyses
sp_by_elev <- data.frame(sp_by_elev)
rownames(sp_by_elev) <- paste(sp_by_elev$full_name, sp_by_elev$elevation, sep = "_")
log_sla <- select(sp_by_elev, log_sla)
log_ldmc <- select(sp_by_elev, log_ldmc)
log_leaf_area <- select(sp_by_elev, log_leaf_area)
log_specific_force <- select(sp_by_elev, log_specific_force) 
stem_density <- select(sp_by_elev, stem_density)
lnc <- select(sp_by_elev, lnc)
log_seed_area <- select(sp_by_elev, log_seed_area)

# Combine into a list for vectorized operations
trait_dfs <- list(log_sla, log_ldmc, log_leaf_area, log_specific_force, stem_density, lnc, log_seed_area)
trait_dfs <- lapply(trait_dfs, scale)
trait_names <- c("log_sla", "log_ldmc", "log_leaf_area", "log_specific_force", "stem_density", "lnc", "log_seed_area")

# Now, apply the distance function to each element of the list
trait_distances <- lapply(trait_dfs, function(x) as.matrix(dist(x)))
names(trait_distances) <- trait_names

# Check what the first few elements in the list look like
# The diagonals should be Zero, pairwise comparisions should be equal in both directions
trait_distances[["log_sla"]][1:5, 1:5]

# Now that we have pariwise distances between all species X site pairs, 
# We can get the degree of turnover between sites

compute_comm_turnover <- function(x, community){
  mean_traitdist_matrix <- matrix(NA, nrow = 5, ncol = 5)
  for (i in 1:5) {
    for (j in 1:5) {
      mean_traitdist_matrix[i,j] <-  mean(x[community[[i]],
                                            community[[j]]], na.rm = T)
    }
  }
  return(as.dist(mean_traitdist_matrix))
}


# Run this function over the trait_dist_list to get a list of distance matrices:
trait_turnovers <- lapply(trait_distances, function(x) 
  compute_comm_turnover(x, species_in_community))
trait_turnovers

# functional dissimilarity with PC Scores ---------
pca <- prcomp(sp_by_elev[,c(6:11)], scale. = T, center = T)
pc_scores <- pca$x

sp_by_elev <- cbind(sp_by_elev, pc_scores)

pca_distance <- as.matrix(dist(data.frame(select(sp_by_elev,13:18))))
pca_turnovers <- compute_comm_turnover(pca_distance, species_in_community)

# Phylo distances ------
phylo_distances <- cophenetic(phylo)


# We need to make a new species in community matrix, because the earlier one has elevation included with species name so that we can get site-specific traits

species_in_community_phylo <- list()
for (i in 1:5) {
  species_in_community_phylo[[i]] <- names(site_by_species[i, site_by_species[i, ] > 0])
}
# Compute Dpw- we can do this with the same compute_comm_turnover() function
# That we had used for functional trait turnovers, since the approach is identical
phylo_diff <- compute_comm_turnover(phylo_distances, species_in_community_phylo)
# Take a look at Dpws
phylo_diff

# We need to also get unifrac distances
# Let's source a unifrac function
source("scripts/functions/gsk_unifrac.R")
phylo_diff_unifrac <- unifrac_gsk(site_by_species, phylo)

# OK, now we have three important sets of difference matrices:
# Taxonomic dissimilarity
tax_diff
trait_turnovers
phylo_diff
phylo_diff_unifrac

# Mantel tests --------
# Correlate Taxon Diff to Elev Diff
tax_mantel <- mantel(tax_diff, elev_diff)

# Correlate Phylo Diff to Elev Diff
phy_mantel <- mantel(phylo_diff, elev_diff)
phy_mantel_unif <- mantel(phylo_diff_unifrac, elev_diff)

# Correlate Trait Diff to Elev Diff
trt_mantel <- lapply(trait_turnovers, function(x) mantel(x, elev_diff))
# Need to update the "xdis" call in the trait_turnover list
for (i in 1:length(trt_mantel)) {
  trt_mantel[[trait_names[i]]]$call$xdis <- trait_names[i]
  print(trt_mantel[[trait_names[i]]]$call$xdis)
}

tax_mantel_temp <- mantel(tax_diff, temp_diff)
phy_mantel_temp <- mantel(phylo_diff, temp_diff)
phy_mantel_unif_temp <- mantel(phylo_diff_unifrac, temp_diff)
trt_mantel_temp <- lapply(trait_turnovers, function(x) mantel(x, temp_diff))
for (i in 1:length(trt_mantel_temp)) {
  trt_mantel_temp[[trait_names[i]]]$call$xdis <- trait_names[i]
  print(trt_mantel_temp[[trait_names[i]]]$call$xdis)
}
pca_mantel <- mantel(pca_turnovers, elev_diff)

# Get all the summary stats into a table
mantel_summary_elevation <- NULL
source("scripts/functions/mantel_summary.R")
mantel_summary_elevation <- t(sapply(trt_mantel, mantel_summary))
mantel_summary_elevation <- rbind(mantel_summary_elevation, mantel_summary(tax_mantel))
mantel_summary_elevation <- rbind(mantel_summary_elevation, mantel_summary(phy_mantel))
mantel_summary_elevation <- rbind(mantel_summary_elevation, mantel_summary(phy_mantel_unif))

# add the temperature mantels
mantel_summary_elevation <- rbind(mantel_summary_elevation, t(sapply(trt_mantel_temp, mantel_summary)))
mantel_summary_elevation <- rbind(mantel_summary_elevation, mantel_summary(tax_mantel_temp))
mantel_summary_elevation <- rbind(mantel_summary_elevation, mantel_summary(phy_mantel_temp))
mantel_summary_elevation <- rbind(mantel_summary_elevation, mantel_summary(phy_mantel_unif_temp))

# Remove rownames since they only exist for traits and are duplicated anyway
rownames(mantel_summary_elevation) <- NULL
View(mantel_summary_elevation)

# Output
if(write_outputs == T) {write.csv(x = mantel_summary_elevation, "tables/mantel_elevation", 
                                  quote = F, row.names = F)}

# Repeat these with Temperature differences, (a proxy of environmental distance)

# Phylogenetic signal in traits ---------
# First we need to get down to a species average across sites...
species_averages <- sp_by_elev %>%  
  group_by(full_name, genus, species, growth_form) %>% 
  summarise(log_sla = mean(log_sla, na.rm = T), 
            log_ldmc = mean(log_ldmc, na.rm = T), 
            log_leaf_area = mean(log_leaf_area, na.rm = T), 
            log_specific_force = mean(log_specific_force, na.rm = T), 
            stem_density = mean(stem_density, na.rm = T),
            lnc = mean(lnc), 
            log_seed_area = mean(log_seed_area, na.rm = T))
rownames(species_averages) <- species_averages$full_name

# Source the function to calculate Blomberg's K using Picante
# This also does a null model randomization
par(default_pars)
source("scripts/functions/k_null_test.R")
# The function is called K_significance
# We can make a list of K values and associated p values
# Let's first do this for all traits except seed projected area.
# We need to do some more munging before working on Seed data because of missing data
K_values <- lapply(trait_names[1:6], 
                   function(x) K_significance(species_averages, x, phylo, 
                                              log = F, reps = reps_for_K, plot = T))
names(K_values) <- trait_names[1:6]

# Now let's get K for seed projected area
# First just make a vector of seed areas to work off of
log_seed_area <- species_averages$log_seed_area
names(log_seed_area) <- rownames(species_averages)

# Figure out which species are missing seed data
taxa_missing_seed <- names(log_seed_area)[which(is.na(log_seed_area))]
# Confirm that we have the right list
taxa_missing_seed

# Drop these from the vector of seed areas
log_seed_area <- log_seed_area[-which(is.na(log_seed_area))]
# Confirm that the length is 101-24 = 77
length(log_seed_area) == 77

# Prune the phylogeny to drop the taxa missing trait data
phylo_for_seed <- drop.tip(phylo, taxa_missing_seed)
length(phylo_for_seed$tip.label) == 77

seed_K <- K_significance(data.frame(log_seed_area), trait = "log_seed_area", 
                         phy = phylo_for_seed, reps = reps_for_K)
# Add this object into the list of Ks we already have
K_values$log_seed_area <- seed_K

# Take a look at all calculated K values
str(K_values)



# Make the phylo+traits figure --------
source("scripts/functions/phy_plottr.R")
# First, reorder the species_averages dataframe. This just helps with plotting
species_averages <- data.frame(species_averages)
species_averages <- species_averages[(phylo$tip.label
                                      [phylo$edge
                                      [phylo$edge[,2]<= Ntip(phylo),2]]),]

# Set some graphing options that help later
pchs <- c(21, 22, 23)
pchvec <- pchs[(species_averages$growth_form)]
colors <- c("black", "grey", "white")
colvec <- colors[(species_averages$growth_form)]
defined_lwd = 0.75
margin_text_size <- 0.65


# Phylo fig ----------------------
if(write_outputs == TRUE) {pdf("figures/phylogeny-fig.pdf", height = 15, width = 21)}
# Set the layout
layout(matrix (c(1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9), 1, 11, byrow = T))
# Plot the phylogeny
par(mar = c(5, 4.1, 4, 1))
plot.phylo(phylo, cex = .90,  label.offset = .03, edge.width = 1.2)
# Add tip labels indicating growth form
tiplabels(pch = pchvec, bg = colvec, cex = 1.25)
# Add a legend for growth forms
legend("topleft", pch = pchs, legend = c(levels(species_averages$growth_form)), 
       bty = "n", horiz = F, cex = 1.2, pt.bg = colors)
n_tips <- length(phylo$tip.label)
# elevation range
par(mar = c(4, 0, 5, 0))
plot(1, type = "n", xlim = c(1, 5), ylim = c(1, length(phylo$tip.label)), 
     bty = "n", axes = F, xlab = "", ylab = "")

mtext(side = 3, line = 1, text = "Distribution", cex = .8)

# Setup
distribution <- data.frame(t((site_by_species > 0)+0))
distribution <- distribution[(phylo$tip.label
                              [phylo$edge
                              [phylo$edge[,2]<= Ntip(phylo),2]]),]
# Make a grey grid for all species
for(ii in 1:n_tips){
  points(x = 1:5, y = rep(ii+1, 5), type = "b", pch = 16, col = "grey", cex = 0.65)
}

# for each species, highlight the site at which they were found
for(ii in 1:n_tips){
  for (jj in 1:5) {
    # print(c(jj, ii))
    points(x = jj, y = ii+1, cex = (distribution[ii, jj]*1), pch = 19)
  }
}

axis(side = 1, at = 1:5, labels = c("30", "500", "800", "2000", "2500"), srt = 45)
axis(side = 3, at = 2:5, labels = c("35", "24", "20", "6"), tick = F, padj = 3)
axis(side = 3, at = 1, labels = c("n species = 45"), tick = F, padj = 3, hadj = .9)
mtext(side = 1, "elevation (m)", line = 2.5, cex = margin_text_size)
box()
# Start plotting traits
phy_plottr(df = species_averages, phy = phylo, trait = "log_sla", col = "grey10", 
           log = F, cex = 1.25, lwd = defined_lwd, 
           title = expression(paste('log'[10],'(SLA)')))
mtext(side = 3, line = 0, text = paste("K =",  round(K_values$log_sla$calc_k, 3)),
      cex = margin_text_size)
mtext(side = 1, expression(paste('cm'^2,'/g')), line = 2.5, cex = margin_text_size)
abline(v = mean(species_averages$log_sla), lwd = 0.7, lty = 3)
box()



# leaf area

phy_plottr(df = species_averages, phy = phylo, trait = "log_leaf_area", col = "grey10",
           log = F, cex = 1.25, lwd = defined_lwd, 
           title = expression(paste('log'[10],'(Leaf Size)')))
mtext(side = 3, line = 0, text = paste("K =",  round(K_values$log_leaf_area$calc_k, 3), "*"),
      cex = margin_text_size)
mtext(side = 1, expression(paste("cm" ^2)), line = 2.5, cex = margin_text_size)
abline(v = mean(species_averages$log_leaf_area), lwd = 0.7, lty = 3)

box()

# Stem Density
phy_plottr(df = species_averages, phy = phylo, trait = "stem_density", col = "grey10",
           log = F, cex = 1.25, lwd = defined_lwd, 
           title = expression(paste('Stem density')))
mtext(side = 3, line = 0, text = paste("K =",  round(K_values$stem_density$calc_k, 3)),
      cex = margin_text_size)
mtext(side = 1, expression(paste('g/','cm'^3)), line = 2.5, cex = margin_text_size)
abline(v = mean(species_averages$stem_density), lwd = 0.7, lty = 3)
box()

# LDMC
phy_plottr(df = species_averages, phy = phylo, trait = "log_ldmc", col = "grey10",
           log = F, cex = 1.25, lwd = defined_lwd, 
           title = expression(paste('log'[10],'(LDMC)')))
mtext(side = 3, line = 0, text = paste("K =",  round(K_values$log_ldmc$calc_k, 3), "*"),
      cex = margin_text_size)
mtext(side = 1, expression(paste('mg/g')), line = 2.5, cex = margin_text_size)
abline(v = mean(species_averages$log_ldmc), lwd = 0.7, lty = 3)
box()

# SFP

phy_plottr(df = species_averages, phy = phylo, trait = "log_specific_force", col = "grey10",
           log = F, cex = 1.25, lwd = defined_lwd, title = expression(paste('log'[10],'(Leaf Toughness)')))
mtext(side = 3, line = 0, text = paste("K =",  round(K_values$log_specific_force$calc_k, 3), "*"),cex = margin_text_size)
mtext(side = 1, expression(paste('N/' ,'m'^2)), line = 2.5, cex = margin_text_size)

abline(v = mean(species_averages$log_specific_force), lwd = 0.7, lty = 3)
box()

# lnc
phy_plottr(df = species_averages, phy = phylo, trait = "lnc", col = "grey10",
           log = F, cex = 1.25, lwd = defined_lwd, title = "Leaf [N]")
mtext(side = 3, line = 0, text = paste("K =",  round(K_values$lnc$calc_k, 3)),
      cex = margin_text_size)
mtext(side = 1, expression(paste('%')), line = 2.5, cex = margin_text_size)
abline(v = mean(species_averages$lnc), lwd = 0.7, lty = 3)
box()

# seed
phy_plottr(df = species_averages, phy = phylo, trait = "log_seed_area", col = "grey10",
           log = F, cex = 1.25, lwd = defined_lwd, 
           title = expression(paste('log'[10],'(Seed Size)')))
mtext(side = 3, line = 0, text = paste("K =",  round(K_values$log_seed_area$calc_k, 3), "*"),
      cex = margin_text_size)
mtext(side = 1, expression(paste(mu, m^2)), line = 2.5, cex = margin_text_size)
abline(v = mean(species_averages$log_seed_area, na.rm = T), lwd = 0.7, lty = 3)
box()


if(write_outputs == TRUE) {dev.off()}


# Regression bootstraps ------
source("scripts/functions/bootstrap_on_regression.R")
source("scripts/functions/aov_bootstrap.R")

traits_boot = c("log_sla", "log_ldmc", "log_leaf_area", "log_specific_force", "lnc", "stem_density", "log_seed_area")

regression_bootstraps = lapply(traits_boot, function(x) regress_bootstrap(df = sp_by_location, trait = x, reps = 10))
aov_bootstraps = lapply(traits_boot, function(x) aov_bootstrap(df = sp_by_location, trait = x))

regression_bootstrap_pvals = sapply(regression_bootstraps, function(x) x$sigpavls)
regressionbootstrap_rs = sapply(regression_bootstraps, function(x) mean(x$rsq))

aov_bootstrap_pvals = sapply(aov_bootstraps, function(x) x$sigpavls)
aov_bootstrap_fs = sapply(aov_bootstraps, function(x) mean(x$fstats))

names(regression_bootstrap_pvals) = traits_boot
names(regressionbootstrap_rs) = traits_boot

names(aov_bootstrap_pvals) = traits_boot
names(aov_bootstrap_fs) = traits_boot

sp_by_location$log_sla
a <- aov(log_sla~as.factor(elevation), data = sp_by_location)
summary(a)
a$effects


## Comparison to Read et al. (2013) results
# elevation to MAT
tmp = data.frame(elevation = c(30,500,800,2000,2500), mat = c(24.6, 22.2, 20.7, 14.5, 11.9))
tmp = merge(tmp,sp_by_location,by="elevation",all=T)

# Read's raw data
rLMA = read.csv("data/READ_lmaraws.csv",as.is=T); head(rLMA) # g/m2
rLNC = read.csv("data/READ_nmraws.csv",as.is=T); head(rLNC) # %

# Filer raw data
rLMA = subset(rLMA,latitude<=24 & elevation<=3000)
rLNC = subset(rLNC,latitude<=24 & elevation<=3000)

plot(lma~elevation,rLMA)
plot(nmass~elevation,rLNC)

# Read's meta-analysis data
Rlma = read.csv("data/READ_lma3.csv",as.is=T); head(Rlma)
Rlnc = read.csv("data/READ_leafn3.csv",as.is=T); head(Rlnc)

# SLA/LMA
cor(1/tmp$sla,tmp$mat) # correlation observed in the CR melastomes (LMA vs. MAT)
# range reported in the meta-analisis: from -0.30 to -0.68

# LNC
cor(tmp$lnc,tmp$mat)
# range reported in the meta-analisis: from -0.?? to 0.??















