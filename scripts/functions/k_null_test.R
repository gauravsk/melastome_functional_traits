# K_significance
# Used to calculate K and assess its significance compared to a null distribution
# Gaurav kandlikar, gkan@umd.edu

K_significance <- function(df, trait, phy, reps=1000, log=TRUE, plot = FALSE){
  
  if(!is.ultrametric(phy)){stop("Phylogeny is not ultrametric!")}
  if(!(trait %in% colnames(df))){stop("Specified trait not present in df!")}
  if(length(setdiff(rownames(df), phy$tip.label) != 0)){stop("rownames of dataframe must match tip labels of phylogeny")}
  
  current_trait <- unlist(df[,trait])
  if (log == TRUE){current_trait <- log10(current_trait)}
  names(current_trait) <- rownames(df)
  
  calculated_k <- Kcalc(x = current_trait, phy = phy)[1,1]
  
  simulated_ks <- numeric(reps+1)  
  
  for (i in 1:reps){
    shuffled_phy <- tipShuffle(phy)
    simulated_ks[i] <- Kcalc(x=current_trait, phy=shuffled_phy) [1,1]
  }
  simulated_ks[reps + 1] = calculated_k

  nbigger <- length(subset(simulated_ks, simulated_ks >= calculated_k))
  prob = nbigger/reps
  
  if (plot == TRUE){
    par(mar=c(3,3,4,3))
    plot(density(simulated_ks), main=paste("Null distribution of K for", trait), xlim=c(0,1.2))
    abline(v = calculated_k, lwd=2, col="darkred")
    mtext(side = 3, text=paste("K = ", round(calculated_k, 3),"; p = ", prob, sep = ""))}
  
  return(list(calc_k = calculated_k, simulated_ks = simulated_ks, p = prob))
}


