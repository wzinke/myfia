#flnm='scan_07_pp_MSE.dat'
#outfl='scan_07_pp_MSE'
#library(mFilter)

# ### TODO !!! This script needs some severe re-work !!!

mse_outlier = function(flnm,outfl,span=0.025)
{
require(KernSmooth)

  dt = read.table(flnm)
  thr_fac = 4
  absMSE = dt$V1
  relMSE = dt$V2
  xval = seq(from=0,by=1,to=length(absMSE)-1)

#   lambdaval=1
#   fltobj=(mFilter(absMSE,filter='HP',drift=F,type='lambda',freq=lambdaval))
#   MSE_trend = fltobj$trend
#   MSE_cycle = fltobj$cycle

  loessfit=loess(absMSE~xval, data.frame(absMSE,xval), span=span, family="symmetric", trace.hat='exact')
  MSE_trend=predict(loessfit,xval)

# repsm = locpoly(xval, absMSE,  bandwidth=span,range.x=range(xval),gridsize=length(xval))
# MSE_trend=repsm$y

  MSE_res = absMSE-MSE_trend
  mad_val = mad(MSE_res)
  med_val = median(MSE_res)
  MSE_cycle = (MSE_res - med_val) / mad_val
  
  uthrv=med_val+thr_fac*mad_val
  lthrv=med_val-thr_fac*mad_val
  thrpos = MSE_res > uthrv | MSE_res < lthrv
  outpos = subset(xval, thrpos)
  num_out=length(outpos)

  rel_uthrv=median(relMSE)+thr_fac*mad(relMSE)
  rel_lthrv=median(relMSE)-thr_fac*mad(relMSE)
  rel_thrpos = relMSE > rel_uthrv | relMSE < rel_lthrv
  reloutpos = subset(xval,rel_thrpos)
  rel_num_out=length(reloutpos)

  write(outpos,    file=paste(outfl,'out.dat',sep='_'), append=F,ncolumns=1)
  write(reloutpos, file=paste(outfl,'relout.dat',sep='_'), append=F,ncolumns=1)
  write(MSE_trend, file=paste(outfl,'trend.dat',sep='_'), append=F,ncolumns=1)
  write(MSE_cycle, file=paste(outfl,'resid.dat',sep='_'), append=F,ncolumns=1)

  write(paste('aMAD: ',mad(absMSE),'\naMED: ',median(absMSE),'\noutlier: ',num_out,'\nresMAD: ',mad_val,'\nresMED: ',med_val,sep=''), file=paste(outfl,'summary.dat',sep='_'), append=F)

  postscript(file=paste(outfl,'.eps',sep=''),paper='a4')
#   postscript(file=paste(outfl,'.eps',sep=''),width = 6*0.39*25, height = 3*0.39*25,paper='special')
  split.screen(c(3,1),erase=T)

  screen(1,new=T)
  par(plt=c(0.05,0.98,0.15,0.95),mgp=c(1.5, 0.4, 0),cex.lab=0.5)
  y_lim = range(c(relMSE))
  plot(xval,relMSE,type='l',ylim=y_lim,xaxs='i',lwd=2,ylab='rel MSE',xlab=paste('MSE Outlier detection: ',flnm,sep=''))
  lines(xval,relMSE,col='blue',lwd=1)
  points(xval[rel_thrpos],relMSE[rel_thrpos],col='green',pch=19)
  abline(v=xval[thrpos],col='red',lty=2)
  abline(h=rel_uthrv,lty=2,col='green')
  abline(h=rel_lthrv,lty=2,col='green')

  screen(2,new=T)
  par(plt=c(0.05,0.98,0.15,0.95),mgp=c(1.5, 0.4, 0),cex.lab=0.5)
  y_lim = range(c(absMSE))
  plot(xval,absMSE,type='l',ylim=y_lim,xaxs='i',lwd=2,ylab='MSE',xlab=paste('MSE Outlier detection: ',flnm,sep=''))
  lines(xval,relMSE,col='blue',lwd=1)
  lines(xval,MSE_trend,col='red',lwd=2)
  points(xval[thrpos],absMSE[thrpos],col='red',pch=19)
  points(xval[rel_thrpos],relMSE[rel_thrpos],col='green',pch=19)
  abline(v=xval[thrpos],col='red',lty=2)
  abline(h=rel_uthrv,lty=2,col='green')
  abline(h=rel_lthrv,lty=2,col='green')

  screen(3,new=T)
  par(plt=c(0.05,0.98,0.15,0.95),mgp=c(1.5, 0.4, 0),cex.lab=0.5)

  plot(xval,MSE_cycle,type='l',xaxs='i',lwd=1.5,xlab='Volume',ylab='normalized MSE residuals [MAD]')
  abline(h=uthrv,lty=2,col='red')
  abline(h=lthrv,lty=2,col='red')
  abline(v=xval[rel_thrpos],col='green',lty=2)
  points(xval[thrpos],MSE_cycle[thrpos],col='red',pch=19)

  text(x=max(xval)/2,y=max(MSE_cycle)-0.02*diff(range(MSE_cycle)),labels=paste('abs Outlier = ',num_out,'    rel Outlier = ',rel_num_out,"    MAD = ",format(mad_val,digits=4),sep=''),font=2,cex=0.8)

  dev.off()

}
