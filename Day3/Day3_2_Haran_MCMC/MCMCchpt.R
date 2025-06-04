## Markov chain Monte Carlo algorithm for a Bayesian (single) change point model
## read in the data
## chptdat = read.table("chpt.dat",header=T)
## Y = chptdat$Ener
## chptdat = read.table("coal.dat",header=T)
## Y = chptdat$Deaths
## Y = read.table("http://personal.psu.edu/muh10/MCMCtut/COUP551_rates.dat", skip = 1)
KGUESS = 10 # our guess for k based on exploratory data analysis
## Note: this function is not written in the most efficient way since its purpose is primarily instructive

## function with default values for length of chain, data set variable, b1, b2
mhsampler = function(NUMIT=1000,dat=Y, b1 = 100, b2 = 100) 
  {
    n = length(dat)
    cat("n=",n,"\n")
    ## set up
    ## NUMIT x 3 matrix to store Markov chain values
    ## each row corresponds to one of 3 parameters in order: theta,lambda,k,b1,b2
    ## each column corresponds to a single state of the Markov chain
    mchain = matrix(NA, 3, NUMIT)
    acc = 0 # count number of accepted proposals (for k only)
    
    ## starting values for Markov chain
    ## This is somewhat arbitrary but any method that produces reasonable values for each parameter is usually adequate.
    ## For instance, we can use approximate prior means or approximate MLEs.
    
    kinit = floor(n/2) # approximately halfway between 1 and n
    mchain[,1] = c(1,1,kinit)
    
    for (i in 2:NUMIT)
      {
        ## most upto date state for each parameter
        currtheta = mchain[1,i-1]
        currlambda = mchain[2,i-1]
        currk = mchain[3,i-1]
        
        ## sample from full conditional distribution of theta (Gibbs update)
        currtheta = rgamma(1,shape=sum(Y[1:currk])+0.5, scale=b1/(currk*b1+1))
        
        ## sample from full conditional distribution of lambda (Gibbs update)
        currlambda = rgamma(1,shape=sum(Y[(currk+1):n])+0.5, scale=b2/((n-currk)*b2+1))
        
        ## sample from full conditional distribution of k (Metropolis-Hastings update)
        propk = sample(x=seq(2,n-1), size=1) # draw one sample at random from uniform{2,..(n-1)}

        ## Metropolis accept-reject step (in log scale)
        logMHratio = sum(Y[1:propk])*log(currtheta)+sum(Y[(propk+1):n])*log(currlambda)-propk*currtheta- (n-propk)*currlambda - (sum(Y[1:currk])*log(currtheta)+sum(Y[(currk+1):n])*log(currlambda)-currk*currtheta- (n-currk)*currlambda)
        
        logalpha = min(0,logMHratio) # alpha = min(1,MHratio)
        if (log(runif(1))<logalpha) # accept if unif(0,1)<alpha, i.e. accept with probability alpha, else stay at current state
          {
            acc = acc + 1 # increment count of accepted proposals
            currk = propk
          }
        
##        currk = KGUESS # if we do not sample k (k fixed) (IF COMMENTED OUT, WE ACCOUNT FOR UNCERTAINTY IN k
        
        ## update chain with new values
        mchain[,i] = c(currtheta,currlambda,currk)
        
      }

    cat("Markov chain algorithm ran for ",NUMIT,"iterations (acc.rate for k=",acc/(NUMIT-1),")\n")
    cat("Parameters are in order: theta, lambda, k\n")
    return(mchain)
  }
