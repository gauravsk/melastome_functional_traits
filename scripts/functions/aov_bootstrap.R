aov_bootstrap <- function(df, trait, reps = 1000) {
  
  elev1 <- unlist(df[df$elevation == 30, trait])
  elev2 <- unlist(df[df$elevation == 500, trait])
  elev3 <- unlist(df[df$elevation == 800, trait])
  elev4 <- unlist(df[df$elevation == 2000, trait])
  elev5 <- unlist(df[df$elevation == 2500, trait])
  
  fstats <- numeric(reps)
  pvals <- numeric(reps)
  for(ii in 1:reps){
    comm1 = sample(x = elev1, size = 6)
    comm2 = sample(x = elev2, size = 6)
    comm3 = sample(x = elev3, size = 6)
    comm4 = sample(x = elev4, size = 6)
    comm5 = sample(x = elev5, size = 6)
    elevs = c(rep(30, 6), rep(500, 6), rep(800, 6), rep(2000, 6), rep(2500, 6))
    df_to_analyse = cbind(elevs, c(comm1, comm2, comm3, comm4, comm5))
    aov_out <- aov(df_to_analyse[,2]~df_to_analyse[,1])
    fstats[ii] <- summary(aov_out)[[1]]$F[1]
    pvals[ii] <- summary(aov_out)[[1]][['Pr(>F)']][1]

  }
  sigpvals = sum(pvals < .05)/length(pvals)
  
  return(list(pvals = pvals, sigpavls = sigpvals, fstats = fstats))
  
}
