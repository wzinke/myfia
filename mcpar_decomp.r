mc.princomp = function(flin, flout=NA, grpnrm=0, nDOF=NA)
{

	dt = read.table(flin)
	attach(dt)
	nvol=length(V1)

	if(is.na(flout)) { flout = flin }
	if(is.na(nDOF))  { nDOF  = length(dt) }
	
# normalize parameter to give them comparable impact
	if(grpnrm == 0) {
		V1 = (V1 - median(V1)) / mad(c(V1))
		V2 = (V2 - median(V2)) / mad(c(V2))
		V3 = (V3 - median(V3)) / mad(c(V3))
		V4 = (V4 - median(V4)) / mad(c(V4))
		V5 = (V5 - median(V5)) / mad(c(V5))
		V6 = (V6 - median(V6)) / mad(c(V6))
	}else{
		V1 = (V1 - median(V1)) / mad(c(V1,V2,V3))
		V2 = (V2 - median(V2)) / mad(c(V1,V2,V3))
		V3 = (V3 - median(V3)) / mad(c(V1,V2,V3))
		V4 = (V4 - median(V4)) / mad(c(V4,V5,V6))
		V5 = (V5 - median(V5)) / mad(c(V4,V5,V6))
		V6 = (V6 - median(V6)) / mad(c(V4,V5,V6))
	}

	if(length(dt) >= 12 & nDOF > 6) {

		if(grpnrm == 0) {
			V7 =  (V7  - median(V7))  / mad(c(V7))
			V8  = (V8  - median(V8))  / mad(c(V8))
			V9  = (V9  - median(V9))  / mad(c(V9))
			V10 = (V10 - median(V10)) / mad(c(V10))
			V11 = (V11 - median(V11)) / mad(c(V11))
			V12 = (V12 - median(V12)) / mad(c(V12))
		}else{
			V7 =  (V7  - median(V7))  / mad(c(V7,V8,V9))
			V8  = (V8  - median(V8))  / mad(c(V7,V8,V9))
			V9  = (V9  - median(V9))  / mad(c(V7,V8,V9))
			V10 = (V10 - median(V10)) / mad(c(V10,V11,V12))
			V11 = (V11 - median(V11)) / mad(c(V10,V11,V12))
			V12 = (V12 - median(V12)) / mad(c(V10,V11,V12))
		}

		dtf=data.frame(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12)   # skip last two entries

	}else{
		dtf=data.frame(V1,V2,V3,V4,V5,V6)
	}


# run PCA
	pdtf=princomp(dtf,cor=F)
	pcobj=pdtf$scores
	
	
# plot components
	postscript(file=paste(flout,'_PC.eps',sep=''),paper='a4')

	par(mfrow=c(2,1),xaxs='i')
	plot(pcobj[,1],type='n',ylim=c(min(pdtf$scores),max(pdtf$scores)),xlim=c(0,nvol),ylab='PC score', xlab='Volume',main='Principal Components of MC parameter')
	lines(pcobj[,1],lwd=0.5,col='black')
	lines(pcobj[,2],lwd=0.5,col='red')
	lines(pcobj[,3],lwd=0.5,col='blue')
	lines(pcobj[,4],lwd=0.5,col='green')
	lines(pcobj[,5],lwd=0.5,col='cyan')
	lines(pcobj[,6],lwd=0.5,col='magenta')

	legend("topleft", c('PC1','PC2','PC3','PC4','PC5','PC6'), col = c('black','red','blue','green','cyan','magenta'),lwd=1.5,ncol = 6,bty='n',cex=0.75)

	plot(pdtf,main='PC variance')

	dev.off()

	par(mfrow=c(1,1))

# write data to file
	write.table(pcobj, file = paste(flout,'_PCmat.dat',sep=''),sep = "    ",row.names=F,col.names=F)

	write.table(pcobj[,1], file = paste(flout,'_PC1.ev',sep=''),sep = " ",row.names=F,col.names=F)
	write.table(pcobj[,2], file = paste(flout,'_PC2.ev',sep=''),sep = " ",row.names=F,col.names=F)
	write.table(pcobj[,3], file = paste(flout,'_PC3.ev',sep=''),sep = " ",row.names=F,col.names=F)
	write.table(pcobj[,4], file = paste(flout,'_PC4.ev',sep=''),sep = " ",row.names=F,col.names=F)
}

