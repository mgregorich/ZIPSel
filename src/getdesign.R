# -----------------------------------------------------------------------------
# Author: GH
# Question: get a permutation of 1:nmax by sample(1:nmax) which fulfills the conditions:
# (1) range of (sample(1:nmax) + nmax:1) == nmax (or nmax-1) (or between nmax/f, nmax*f)
# (2) all values of (sample(1:nmax) + nmax:1) are unique
# -----------------------------------------------------------------------------

rm(list=ls())
nmax <- 25   ## number of peptides
f <- 1.5    ## tolerance factor (can be 1 for small values of nmax, should be >1 for larger values)

set.seed(7123981)
nsim <- 100000 # can be smaller for small nmax

b <- matrix(0,nsim,nmax)
good <- (1:nsim)*0
luni <- diffrange <- (1:nsim)*0

for(i in 1:nsim){
  a <-c(seq(nmax,1, by=-2),seq(1,nmax, by=2))  # these are the ranks of beta_u for peptides 1:nmax
  b[i,] <- sample(1:nmax) # these are the ranks of beta_d for peptides 1:nmax
  c1 <- a+b[i,] # ranksum(beta_u, beta_d) for a particular permutation
  diffrange[i] <- diff(range(c1))  # compute range of ranksum
  luni[i] <- length(unique(c1))  # compute number of different ranksums
}
corr <- 0 # correction for odd nmax (number of peptides)
if(nmax/2 != floor(nmax/2)) corr <- 1 # if odd then corr = 1

# (1) first condition: range of ranksum should be within [nmax/f, nmax*f], where f is a tolerance factor
f1 <- 2
cond1 <- (diffrange <= (nmax-corr)*f1) & (diffrange >= (nmax-corr)/f1)
sum(cond1)

# (2) second condition: all ranksums must be unique (no equal ranksums), harder to achieve! can be turned off
# good <- good * (luni==nmax)
f2 <- 1.25
cond2 <- (luni>=nmax/f2) & luni==max(luni)
sum(cond2)

good <- cond1 * cond2
sum(good)

bgood <- b[good==1,,drop=FALSE]
bgood <- bgood[order(-abs(diffrange[good]-(nmax-corr))),,drop=FALSE]
# 'good' permutations; better ones (in terms of range is closer to nmax) on top

# transform to character and show only unique permutations (unique is easier with character)
# these would be the permutations of the rankings for beta_d (assuming the ranking for beta_u is nmax:1)
chgood <- unique(apply(bgood,1,function(X) paste(X, collapse="-")))
chgood

luni[good]
diffrange[good]
