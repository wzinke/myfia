## infl = "/home/zinke/Dropbox/BRAVO/smoothed/ATTENTION/bca_ab_p_BRAVO_multimediation.nii.gz"
## ofl  = "/home/zinke/Dropbox/BRAVO/smoothed/thr_slices/ATTENTION/bca_ab_p_tst"

vol_MCC = function(infl, ofl=NULL, method='holm', mask=NULL){
# apply multiple comparison correction on a 3D data file (nifti format).
# currently implemented to just utilize p.adjust

if(is.null(ofl)){  ofl = paste(unlist(strsplit(infl,"[.]"))[1],'_mcc.nii.gz',sep='')  }

# # the fslr package might be handy as well
# require('fslr')
# # set fslr options
# options(fsl.path="/usr/local/fsl")
# have.fsl() # check
# options(fsl.outputtype = "NIFTI_GZ")

# need the oro.nifti package to import nifti files into R
require('oro.nifti')

# read in the data file
img = readNIfTI(infl)

img[img==0] = NA # set voxels outside of mask to NA - correct this before saving to avoid FSL trouble

if(!is.null(mask)){
  maskfl = readNIfTI(mask)
  img(maskfl==0) = NA
}

img2 = p.adjust(img, method=method) # this returns a vector

# transform to 3D volume again (check for possible dimension flips)
oimg = array(img2, dim=dim(img))

# convert to nifti format
oimg2 = nifti(oimg, value=img)

writeNIfTI(oimg2, ofl, onefile=TRUE, gzipped=TRUE)

}