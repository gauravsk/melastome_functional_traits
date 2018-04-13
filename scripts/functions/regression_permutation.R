regression_permutation <- function(df, trait, reps = 1000) {
  
  trait <- unlist(df[,trait])
  elev <- unlist(df[,"elevation"])
  
  rsq <- numeric(reps)
  pvals <- numeric(reps)
  for(ii in 1:reps){
    df_to_analyse <- cbind(sample(elev), trait)
    net_lm <- lm(df_to_analyse[,2]~df_to_analyse[,1])
    rsq[ii] <- (summary(net_lm)$r.squared)
    pvals[ii] <- anova(net_lm)$'Pr(>F)'[1]
  }
  sigpvals = sum(pvals < .05)/length(pvals)
  true_lm <- lm(trait~elev)
  true_rsq <- (summary(true_lm)$r.squared)
  pperm <- 1-(sum(true_rsq > rsq))/reps
  
  return(list(pvals = pvals, rsq = rsq, sigpavls = sigpvals, true_rsq = true_rsq, pperm = pperm))
  
}

