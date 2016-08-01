# Hand-coded calculation of Binary Jaccard Distance
# Recall that Jaccard is computed as
# 1 - (Number of Shared species)/(Number unique to Comm1 + Number unique to Comm2 + Number Shared)
# site_by_species needs to be a community matrix.
# A distance matrix is returned. 
# The output of this function should be identical to vegan::vegdist(method = "jaccard", binary = T)
jaccard <- function(site_by_species){
  species_by_community <- list()
  
  # Make a list of species in each community
  for (i in 1:nrow(site_by_species)) {
    species_by_community[[i]] <- paste(names(site_by_species[i, site_by_species[i, ] > 0]))
  }
  
  distance_matrix <- matrix(nrow  = length(species_by_community), ncol = length(species_by_community))
  rownames(distance_matrix) <- rownames(site_by_species)
  colnames(distance_matrix) <- rownames(site_by_species)

  for (i in 1:nrow(distance_matrix)){
    for (j in 1:ncol(distance_matrix)) {
      in_both     <- intersect(species_by_community[[i]], species_by_community[[j]])
      in_first    <- sum(!(species_by_community[[i]] %in% species_by_community[[j]]))
      in_second   <- sum(!(species_by_community[[j]] %in% species_by_community[[i]]))
      
      # Use those values to get the distance
      jac_dif     <- 1-length(in_both)/(in_first+in_second+length(in_both))
      
      # Put that distance into the returning matrix
      distance_matrix[i,j] <- jac_dif
    }
  }
  return(as.dist(distance_matrix))
}