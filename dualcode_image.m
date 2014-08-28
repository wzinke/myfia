function oimg = dualcode_image(bckimg, statmap, betamap, sigmap, betarng)
% Create overlay images that display the effect size (parameter estimate)
% as color map, and use the image alpha value as hue to represent statistical
% significance. A binary map selects significant voxels that are highlighted
% as contour and displayed as solid color.

% based on the dualcodeImage function of E.A. Allen (http://mialab.mrn.org/datavis/)

% ToDo: Use nifti files as input instead of image files.

if(~exist('bckimg','file'))
    error('No valid background image specified!');
else
    BCK = imread(bckimg);
end

if(~exist('statmap','file'))
    error('No valid statistical image specified!');
else
    ZMAP = imread(statmap);
end

if(~exist('betamap','file'))
    error('No valid effect size image specified!');
else
    BMAP = imread(betamap);
end

if(~exist('sigmap','file'))
    error('No valid threshold image specified!');
else
    SMAP = imread(sigmap);
end


% Set the Min/Max values for hue coding
% Hue codes the effect size
if(~exist('betarng','var') || isempty(betarng))
    absmax = max(abs(BMAP(:)));
    H_range = [-absmax absmax]; % The colormap is symmetric around zero
elseif(length(betarng) == 1)
    absmax = abs(betarng);
    H_range = [-absmax absmax];
end

