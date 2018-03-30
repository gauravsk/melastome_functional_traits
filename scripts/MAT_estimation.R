### Estimation of mean annual temperature for each site (elevation) based on 
# data reported by Clark et al. (2015)

## original data
elev = c(129,188,417,579,954,1448,1994,2397,2829) # elevation in meters
mat = c(25,23.1,22.3,21.6,20,17.8,14.5,12.1,10.4) # mean temperature in °C

mod = lm(mat~elev) # linear model
mod

plot(mat~elev)
abline(mod,col=2)


## sites sampled in our study
obs = c(30,500,800,2000,2500) # sampled elevations
est = mod$coef[1]+mod$coef[2]*obs # estimated MAT
round(est,1)
diff(range(est)) # MAT span across the transect
