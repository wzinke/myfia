
# flnm='RvsL_brainacc.dat'
# flnm='SvsC_brainacc.dat'
AccHist = function(flnm, onm, chlvl=0.5){
   library(mixtools)

    acc = as.numeric(scan(flnm))

    postscript(file=onm,width = 0.39*20, height = 0.39*20,paper='special',horizontal=F)

    # mixmdl = spEMsymloc(acc, c(0.5,0.6), maxiter = 1000)

    par(xaxs="i",las=1)

#     plot(density(acc, na.rm=T, bw="SJ"), type='n',
#     xlab="classification accuracy", xlim=c(0,1), main=flnm)

    hist(acc, breaks=seq(from=-0.0125,to=1.0125,by=0.025), add=F, prob=TRUE, col="black", border="white",  xlab="classification accuracy", xlim=c(0,1), main=flnm)

    lines(density(acc,na.rm=T),col="blue",lwd=4)
    abline(v=0.5,lty=2,col="orange",lwd=4)
    abline(v=median(acc),lty=2,col="blue",lwd=4)

    tryCatch({
#         mixmdl = normalmixEM(acc, mu=c(0.6,chlvl), maxit=10000, mean.constr=c(chlvl,">=0.5b"), maxrestarts=100, lambda=c(0.3,0.7))
        mixmdl = normalmixEM(acc, mu=c(chlvl,1.25*chlvl), maxit=10000, mean.constr=c(chlvl,NA), fast=FALSE, maxrestarts=100, lambda=c(0.3,0.7))

        amp1   = mixmdl$lambda[1]
        mu1    = mixmdl$mu[1]
        sigma1 = mixmdl$sigma[1]
        curve(amp1*dnorm(x,mu1,sigma1),add=TRUE,lty=1,col="red",lwd=5)

        amp2   = mixmdl$lambda[2]
        mu2    = mixmdl$mu[2]
        sigma2 = mixmdl$sigma[2]
        curve(amp2*dnorm(x,mu2,sigma2),add=TRUE,lty=1,col="green",lwd=5)

        abline(v=mu1,lty=2,col="red",  lwd=4)
        abline(v=mu2,lty=2,col="green",lwd=4)
    })
    # abline(v=0.4,lty=2,col="green",lwd=2.5)
    # abline(v=0.6,lty=2,col="green",lwd=2.5)

    dev.off()
}
