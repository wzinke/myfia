function oimg = dualcode_image(bckimg, statmap, betamap, sigmap)
% Create overlay images that display the effect size (parameter estimate)
% as color map, and use the image alpha value as hue to represent statistical
% significance. A binary map selects significant voxels that are highlighted
% as contour and displayed as solid color.

% based on the dualcodeImage function of E.A. Allen (http://mialab.mrn.org/datavis/)