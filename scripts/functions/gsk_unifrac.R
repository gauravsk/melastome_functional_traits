
unifrac_gsk <- function(site_by_species, phylo){
  # Set up the distance matrix to be returned
  distance_matrix <- matrix(nrow  = nrow(site_by_species), ncol = nrow(site_by_species))
  rownames(distance_matrix) <- rownames(site_by_species)
  colnames(distance_matrix) <- rownames(site_by_species)

  for (i in 1:nrow(distance_matrix)){
    for (j in 1:ncol(distance_matrix)) {
      # Total Phylogenetic Diversity in community 1
      pd_first <- pd(site_by_species, phylo)[i, 1]
      # Total Phylogenetic diversity in community 2
      pd_second <- pd(site_by_species, phylo)[j, 1]
      # Total phylogenetic diversity across both communities
      pd_both <- pd(t(as.matrix(site_by_species[i,])) + 
                      t(as.matrix(site_by_species[j,])), phylo)[1,1]
      # The portion of the total phylo diversity that is unique to one of the two communities
      pd_shared <- (pd_first + pd_second) - pd_both
      # Unifrac is calculated as the 
      # (phylo div unique to one of the two)/(total phylo div)
      u_frac <- (pd_both-pd_shared)/pd_both
      distance_matrix[i,j] <- u_frac

      # Put that distance into the returning matrix
      distance_matrix[i,j] <- u_frac
    }
  }
  return(as.dist(distance_matrix))
}