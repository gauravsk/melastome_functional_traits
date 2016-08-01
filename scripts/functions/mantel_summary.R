# Function to summarize the results of a Mantel correlation
# x needs to be an object of the class "Mantel"
mantel_summary <- function(x){
  if(!class(x) == "mantel") {stop("This function's input needs to be of class 'mantel'")}

  to_return <- as.character(c("call" = x$call$xdis, "ydis" = x$call$ydis, "Mantel R" = x$statistic, "p" = x$signif))
  names(to_return) <- c("xdis", "ydis", "Mantel R", "p")
  return(to_return)
}