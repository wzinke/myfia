# This function reads in the motion parameter file created by
# motion correction routines.
# based on the file format a distinction is made between 3dvolreg and mcflirt.
# Plots are saved as eps-file.
#
#
# wolf zinke, Nov/Dec 2007

mc.plots = function(flnm,epsfl)
{
    mcpar = read.table(flnm,header=F)

    if (length(mcpar) ==  6 | length(mcpar) ==  14 ){
#         xrot = mcpar$V1 * (180/pi)
#         yrot = mcpar$V2 * (180/pi)
#         zrot = mcpar$V3 * (180/pi)
        xrot = mcpar$V1
        yrot = mcpar$V2
        zrot = mcpar$V3
        xtra = mcpar$V4
        ytra = mcpar$V5
        ztra = mcpar$V6
    }
    else{
#         xrot = mcpar$V1 * (180/pi)
#         yrot = mcpar$V2 * (180/pi)
#         zrot = mcpar$V3 * (180/pi)
        xrot = mcpar$V1
        yrot = mcpar$V2
        zrot = mcpar$V3
        xtra = mcpar$V5
        ytra = mcpar$V6
        ztra = mcpar$V7
    }


# get summary for translations and rotations

med_traX = median(xtra)
mad_traX =    mad(xtra)
med_traY = median(ytra)
mad_traY =    mad(ytra)
med_traZ = median(ztra)
mad_traZ =    mad(ztra)

med_rotX = median(xrot)
mad_rotX =    mad(xrot)
med_rotY = median(yrot)
mad_rotY =    mad(yrot)
med_rotZ = median(zrot)
mad_rotZ =    mad(zrot)


mn_traX = mean(xtra)
sd_traX =   sd(xtra)
mn_traY = mean(ytra)
sd_traY =   sd(ytra)
mn_traZ = mean(ztra)
sd_traZ =   sd(ztra)

mn_rotX = mean(xrot)
sd_rotX =   sd(xrot)
mn_rotY = mean(yrot)
sd_rotY =   sd(yrot)
mn_rotZ = mean(zrot)
sd_rotZ =   sd(zrot)

write(paste('MEDIAN_X_tra: ',med_traX,'\nMAD_X_tra: ',mad_traX,'\nMEDIAN_Y_tra: ',med_traY,'\nMAD_Y_tra: ',mad_traY,'\nMEDIAN_Z_tra: ',med_traZ,'\nMAD_Z_tra: ',mad_traZ, '\nMEDIAN_X_rot: ',med_rotX,'\nMAD_X_rot: ',mad_rotX,'\nMEDIAN_Y_rot: ',med_rotY,'\nMAD_Y_rot: ',mad_rotY,'\nMEDIAN_Z_rot: ',med_rotZ,'\nMAD_Z_rot: ',mad_rotZ,   '\n##########\nMEAN_X_tra: ',mn_traX,'\nSD_X_tra: ',sd_traX,'\nMEAN_Y_tra: ',mn_traY,'\nSD_Y_tra: ',sd_traY,'\nMEAN_Z_tra: ',mn_traZ,'\nSD_Z_tra: ',sd_traZ, '\nMEAN_X_rot: ',mn_rotX,'\nSD_X_rot: ',sd_rotX,'\nMEAN_Y_rot: ',mn_rotY,'\nSD_Y_rot: ',sd_rotY,'\nMEAN_Z_rot: ',mn_rotZ,'\nSD_Z_rot: ',sd_rotZ,sep=''), file=paste(epsfl,'summary.dat',sep='_'), append=F)

# make plots
    max_rot = max(c(zrot,yrot,xrot))
    min_rot = min(c(zrot,yrot,xrot))
    max_tra = max(c(ztra,ytra,xtra))
    min_tra = min(c(ztra,ytra,xtra))

    nvol = length(xtra)
    xval = seq(from=1,to=nvol)

	postscript(file=paste(epsfl,'_rot_tra.eps',sep=''),paper='a4')

	par(mfrow=c(2,1),xaxs='i')
	plot(xval,zrot,ylim=c(min_rot,max_rot),type='l',xlim=c(1,nvol),col='blue',lwd=1.5,main='Rotations',ylab='Rotation [deg]',xlab='Volume #',cex.lab=1.5,cex.main=2)
	lines(xval,xrot,col='red',lwd=1.5)
	lines(xval,yrot,col='darkgreen',lwd=1.5)
	abline(0,0,col='black')
	legend("topleft", c('X-rot','Y-rot','Z-rot'), col = c('red','darkgreen','blue'),lwd=2,lty = c(1, 1, 1),inset=0.01,ncol = 3,bty='n')

	plot(xval,ztra,ylim=c(min_tra,max_tra),type='l',xlim=c(1,nvol),col='blue',lwd=1.5,main='Translations',ylab='Translation [mm]',xlab='Volume #',cex.lab=1.5,cex.main=2)
	lines(xval,xtra,col='red',lwd=1.5)
	lines(xval,ytra,col='darkgreen',lwd=1.5)
	abline(0,0,col='black')
	legend("topleft", c('X-tra','Y-tra','Z-tra'), col = c('red','darkgreen','blue'),lwd=2,lty = c(1, 1, 1),inset=0.01,ncol = 3,bty='n')

	dev.off()

	par(mfrow=c(1,1))

    if(length(mcpar) ==  14){

        xsca = mcpar$V7
        ysca = mcpar$V8
        zsca = mcpar$V9
        xyskew = mcpar$V10
        xzskew = mcpar$V11
        yzskew = mcpar$V12

        max_sca = max(c(zsca,ysca,xsca))
        min_sca = min(c(zsca,ysca,xsca))
        max_skew = max(c(xyskew,xzskew,yzskew))
        min_skew = min(c(xyskew,xzskew,yzskew))

        postscript(file=paste(epsfl,'_sca_skew.eps',sep=''),paper='a4')

        par(mfrow=c(2,1))

        plot(xval,zsca,ylim=c(min_sca,max_sca),type='l',xlim=c(1,nvol),col='blue',lwd=1.5,main='Scalings',ylab='Scaling',xlab='Volume #',cex.lab=1.5,cex.main=2)
        lines(xval,xsca,col='red',lwd=1.5)
        lines(xval,ysca,col='darkgreen',lwd=1.5)
        abline(0,0,col='black')
        legend("topleft", c('X-scale','Y-scale','Z-scale'), col = c('red','darkgreen','blue'),lwd=2,lty = c(1, 1, 1),inset=0.01,ncol = 3,bty='n')

        plot(xval,yzskew,ylim=c(min_skew,max_skew),type='l',xlim=c(1,nvol),col='blue',lwd=1.5,main='Skews',ylab='Skew',xlab='Volume #',cex.lab=1.5,cex.main=2)
        lines(xval,xyskew,col='red',lwd=1.5)
        lines(xval,xzskew,col='darkgreen',lwd=1.5)
        abline(0,0,col='black')
        legend("topleft", c('XY-skew','XZ-skew','YZ-skew'), col = c('red','darkgreen','blue'),lwd=2,lty = c(1, 1, 1),inset=0.01,ncol = 3,bty='n')

    dev.off()

    par(mfrow=c(1,1))

}

}

