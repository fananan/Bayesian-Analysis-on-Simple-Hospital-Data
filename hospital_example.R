source("DBDA2E-utilities.R") 
require(rjags)              

fileNameRoot="Jags-ExampleScript" # For output file names.

# Load the data:
data = data.frame(
  n = c(27,148,119,810,211,196,148,215,207,97,256,360),
  y = c(0,18,8,46,8,13,9,31,14,8,29,24)
) 
y = data$y       
n = data$n
Ntotal = length(y)  # Compute the total number of flips.
dataList = list(    # Put the information into a list.
  y = y ,
  Ntotal = Ntotal,
  n = n
)

#beta distributions
omega = .4
mu = .5
kappa = 50
hist(rbeta(10000,omega*(kappa-2)+1,(1-omega)*(kappa-2)+1))
hist(rbeta(10000,mu*kappa,(1-mu)*kappa))

# Define the model:
modelString = "
model {
for ( i in 1:Ntotal ) {
  y[i] ~ dbin(theta[i],n[i]) #you rarely can challenge this part, dead/live, binomial
  theta[i] ~ dbeta(omega*(kappa-2)+1,(1-omega)*(kappa-2)+1)  #we make the subjective knowledge here, we choose beta distribution
}
omega ~ dbeta(1,1) #uniform vague
kappa <- kappaMinusTwo + 2
kappaMinusTwo ~ dgamma(.01,.01)
}
"
writeLines( modelString , con="TEMPmodel.txt" )

#the other model
model{
  for (i in 1:Ntotal){
    y[i]~dbin(theta[i],n[i])
    logit(theta[i])<-b[i] #trick: you spread theta out to infinity, then you are not limited only in beta, since the range is changed; log(theta/(1-theta))
    b[i]~dnorm(mu,tau)
  }
  mu~dnorm(0)
  sigma<-1/sqrt(tau)
  tau~dgamma(0,10)
}


modelString = "
model {
for ( i in 1:Ntotal ) {
  y[i] ~ dbin(theta[i],n[i])
  theta[i] ~ dbeta(omega*kappa,(1-omega)*kappa)
}
omega ~ dbeta(1,1) #uniform vague
kappa ~ dgamma(.01,.01)
}
"# close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Initialize the chains based on MLE of data.
# Option: Use single initial value for all chains:
#  thetaInit = sum(y)/length(y)
#  initsList = list( theta=thetaInit )
# Option: Use function that generates random values for each chain:
initsList = function() {
  return(list(omega=runif(1),kappaMinusTwo=runif(1),theta =runif(Ntotal)))
}

# Run the chains:
jagsModel = jags.model( file="TEMPmodel.txt" , data=dataList , inits=initsList , 
                        n.chains=3 , n.adapt=500 )
update( jagsModel , n.iter=500 )
codaSamples = coda.samples( jagsModel , variable.names=c("omega","kappa","theta") ,
                            n.iter=3334 )

# Examine the chains:
# Convergence diagnostics:
#diagMCMC( codaObject=codaSamples , parName="theta" )
#saveGraph( file=paste0(fileNameRoot,"ThetaDiag") , type="eps" )
plot(codaSamples)
summary(codaSamples)
# Posterior descriptives:
par(mfrow=c(1,1))
samps = do.call(rbind,codaSamples)
hist(samps[,'omega'])
par(mfrow=c(2,1))
hist(samps[,'theta[1]'],xlim=c(0,.2))
hist(samps[,'theta[5]'],xlim=c(0,.2))
plot(apply(samps[,3:14],2,mean),main='model estimates',ylim=c(0,.2),ylab='')
plot(y/n,main='raw estimates',ylim=c(0,.2),ylab='')
par(mfrow=c(1,1))
hist(samps[,'theta[1]']-samps[,'theta[5]'],main='A-E')
hist(samps[,'theta[8]']-samps[,'theta[5]'],main='H-E')
hist(samps[,'theta[8]']-samps[,'theta[1]'],main='H-A')
mean(samps[,'omega'])
sum(y)/sum(n)

plotPost( codaSamples[,"y_new"] , main="Posterior Predictive Distribution" , xlab=bquote(y_new) )
# Re-plot with different annotations:
plotPost( codaSamples[,"diff"] , main="mu1-mu2" , xlab=bquote(diff) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )
