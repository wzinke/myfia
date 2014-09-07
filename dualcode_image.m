function dualcode_image(bckimg, statmap, betamap, sigmap, slcs, sldim, betarng, alpharng, bckrng, ofl, ip, scl)
% Create overlay images that display the effect size (parameter estimate)
% as color map, and use the image alpha value as hue to represent statistical
% significance. A binary map selects significant voxels that are highlighted
% as contour and displayed as solid color.

% based on the dualcodeImage function of E. A. Allen (http://mialab.mrn.org/datavis/)
%     Elena A. Allen, Erik B. Erhardt, & Vince D. Calhoun (2012)
%     Data Visualization in the Neurosciences: Overcoming the Curse of Dimensionality
%     Neuron 74, 603 - 608
%
% wolf zinke, Sep. 2014

%% check inputs and load nifti files
%% get slices
if(~exist('ofl','var'))
    ofl = [];
end

% scale factor [obsolet: interpolation steps (and thus increase image size by this factor)}
if(~exist('ip','var') || isempty(ip))
    ip = 1;
end

% image scale (increase size without interpolation)
if(~exist('scl','var') || isempty(scl))
    if(ip==0)
        scl = 2;
    else
        scl = 1;
    end
end

scrres = get(0,'ScreenSize');

XYctr = round(scrres([3,4])./2);

if(isnumeric(bckimg))
    BCK = bckimg;
elseif(exist(bckimg,'file'))
    chk_FSL;
    BCK = read_avw(bckimg);
else
    error('No valid background image specified!');
end

dims = size(BCK);

if(isnumeric(bckimg))
    ZMAP = statmap;
elseif(exist(statmap,'file'))
    chk_FSL;
    ZMAP = read_avw(statmap);
else
    error('No valid statistical image specified!');
end

if(isnumeric(bckimg))
    BMAP = betamap;
elseif(exist(betamap,'file'))
    chk_FSL;
    BMAP = read_avw(betamap);
else
    error('No valid effect size image specified!');
end

if(~exist('sigmap','var') || isempty(sigmap))
    SMAP = abs(ZMAP) > 2.3;
elseif(isnumeric(sigmap) && length(sigmap) == 1)
    SMAP = abs(ZMAP) > sigmap;
elseif(isnumeric(sigmap))
    SMAP = sigmap;
elseif(exist(sigmap,'file'))
     SMAP = read_avw(sigmap);
else
    error('No valid threshold image specified!');
end

%% get slices
if(~exist('sldim','var') || isempty(sldim))
    sldim = 'z';
end

if(~exist('slcs','var') || isempty(slcs))
    slcs = 0.5;
end

switch lower(sldim)
    case 'x'
        if(any(slcs < 1))
            slcs = ceil(dims(1) * slcs);
        end
        x = slcs;
        if(any(slcs > dims(1)))
            error('Slice exceeds number of slices in volume');
        end
        
        y=1:dims(2);
        z=1:dims(3);
    case 'y'
        if(any(slcs < 1))
            slcs = ceil(dims(2) * slcs);
        end
        y = slcs;
        if(any(slcs > dims(2)))
            error('Slice exceeds number of slices in volume');
        end
        x=1:dims(1);
        z=1:dims(3);
    case 'z'
        if(any(slcs < 1))
            slcs = ceil(dims(3) * slcs);
        end
        z = slcs;
            if(any(slcs > dims(3)))
                error('Slice exceeds number of slices in volume');
            end
        x=1:dims(1);
        y=1:dims(2);
    otherwise
        error('wrong slice dimension specified!');
end

%% set plot ranges
% Set the Min/Max T-values for alpha coding
% opaqueness codes statistical t/z value
if(~exist('bckrng','var') || isempty(bckrng))
    bckvals = BCK(x,y,z);
    bckvals = bckvals(:);
    bckvals(bckvals == 0) = [];
    B_range  = prctile(bckvals,[2 98]);
elseif(length(bckrng) == 1)
    B_range = [0 bckrng];
else
    B_range = bckrng;
end

% Set the Min/Max values for hue coding
% Hue codes the effect size
if(~exist('betarng','var') || isempty(betarng))
    bvals = BMAP(x,y,z);
    bvals = abs(bvals(:));
    absmax  = max(prctile(bvals(bvals>0),[1, 99]));
%     absmax  = max(abs(bvals(:)));
    H_range = [-absmax absmax]; % The colormap is symmetric around zero
elseif(length(betarng) == 1)
    absmax  = abs(betarng);
    H_range = [-absmax absmax];
else
    H_range = betarng;
end

% Set the Min/Max T-values for alpha coding
% opaqueness codes statistical t/z value
if(~exist('alpharng','var') || isempty(alpharng))
    alphvals = abs(ZMAP(x,y,z));
    alpharng  = max(prctile(alphvals(:),[99]));
    A_range = [0 alpharng]; % The colormap is symmetric around zero
elseif(length(alpharng) == 1)
    A_range = [0 alpharng];
else
    A_range = alpharng;
end

for(i = 1:length(slcs))
    eval([lower(sldim), ' = ',int2str(slcs(i)),';']); % careful, this redefines one dimension vector!
    %% Transform the underlay and beta map to RGB values, based on specified colormaps
    if(ip > 1)
        % ipBCK    = interp2(rot90(squeeze(BCK(x,y,z))),  ip,'cubic');
        % ipPE     = interp2(rot90(squeeze(BMAP(x,y,z))), ip,'cubic');
        % alphamap = interp2(rot90(squeeze(ZMAP(x,y,z))), ip,'cubic');
        % ipSIG    = interp2(rot90(squeeze(SMAP(x,y,z))), ip,'cubic');

        ipBCK    = imresize(rot90(squeeze(BCK(x,y,z))),  ip, 'bicubic');
        ipPE     = imresize(rot90(squeeze(BMAP(x,y,z))), ip, 'bicubic');
        alphamap = imresize(abs(rot90(squeeze(ZMAP(x,y,z)))), ip, 'bicubic');
        ipSIG    = imresize(rot90(squeeze(SMAP(x,y,z))), ip, 'bicubic');

        ipSIG(ipSIG > 0.5) = 1; ipSIG(ipSIG <= 0.5) = 0;
    else
        ipBCK    = rot90(squeeze(BCK( x,y,z)));
        ipPE     = rot90(squeeze(BMAP(x,y,z)));
        alphamap = abs(rot90(squeeze(ZMAP(x,y,z))));
        ipSIG    = rot90(squeeze(SMAP(x,y,z)));
    end

    CM_under = bone(256);      % colormap for the underlay (anatomical)
    CM_over  = twowaycol(256); % color map for the effect size

    U_RGB = convert_to_RGB(ipBCK, CM_under, B_range);
    O_RGB = convert_to_RGB(ipPE,  CM_over,  H_range);

    % Use the T/Z-statistics to create an alpha map (which must be in [0,1])
    alphamap(alphamap > A_range(2)) = A_range(2);
    alphamap(alphamap < A_range(1)) = 0;
    alphamap = alphamap/A_range(2);

    %% plot data
    % Make a figure and set of axes
    imsz = size(ipBCK) .* scl;
    figpos = [ XYctr(1) - ceil(imsz(1)/2),  XYctr(2) - ceil(imsz(2)/2), imsz(1), imsz(2)];
    F = figure('Color', 'k', 'Units', 'pixels', 'Position', figpos);
    axes('Position', [0 0 1 1]);

    % Plot the underlay
    image(U_RGB);
    hold on;

    % Now, add the Beta difference map as an overlay
    layer2 = imagesc(O_RGB);

    % Adjust the alpha values of the overlay
    set(layer2, 'alphaData', alphamap);
    alpha(layer2,alphamap);

    % Add some (black) contours to annotate nominal significance
    if(sum(ipSIG(:)) > 0)
        contour(ipSIG, 1, 'k', 'LineWidth', (ip+1)/2);
    end
    
    axis off;
    axis image;

    %% save file
    if(~isempty(ofl))
        cfl = [ofl,'_', upper(sldim),int2str(slcs(i)),'.png'];
        hgexport(F, cfl, hgexport('factorystyle'), 'Format', 'png');
        % print(F, cfl,'-dpng','-r0');
        crop(cfl, 0, 0); % get rid of the white margins Matlab added to the image frame
        close(F);
    end
end

%% 4. Create a 2D colorbar for the dual-coded overlay
%--------------------------------------------------------------------------
figpos = [ XYctr(1) - ceil(0.85*imsz(1)),  XYctr(2) - ceil(imsz(2)/2), round(0.35*imsz(1)), imsz(2)];

G = figure('color', 'k', 'Position', figpos);

x = linspace(A_range(1), A_range(2), 256); % range in alpha (abs(t/z-stats))
y = linspace(H_range(1), H_range(2), size(CM_over,1)); % range in hue (beta weight)
[X,Y] = meshgrid(x,y); % Transform into a 2D matrix

imagesc(x,y,Y);
axis xy; % Plot the colorbar

colormap(CM_over);
alpha(X);
alpha('scaled');

set(gca, 'Xcolor', 'w', 'Ycolor', 'w', 'FontSize', 10, 'LineWidth',2)
set(gca, 'YAxisLocation', 'right', 'TickDir', 'out')
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'Color',[ 0.65 0.65 0.65])
% set(gca,'Color',[1 1 1])

axis tight

xlabel('z value');
ylabel('parameter estimate');

if(~isempty(ofl))
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'PaperPositionMode','auto')
    cfl = [ofl,'_', upper(sldim),'_colbar.png'];
    print(G,cfl,'-dpng','-r0');
    % cfl = [ofl,'_', upper(sldim),'_colbar.eps'];
    % print(G,cfl,'-depsc2','-r0');
    % hgexport(G, cfl, hgexport('factorystyle'), 'Format', 'png');
    close(G);
end

%--------------------------------------------------------------------------
%% check for matlab installation
function chk_FSL
    if(~exist('read_avw','file'))
        [~, fslpth] = system('echo $FSLDIR');

        fslmat = fullfile(fslpth(1:end-1), 'etc/matlab');
        if(~exist(fslmat,'dir'))
            error('FSL ist not correctly installed on the system!')
        else
            addpath(fslmat);
        end
    end

%% Helper function: convert_to_RGB
% this function was provided by E. A. Allen with the dualcodeimage code
% (http://mialab.mrn.org/datavis/)
function IMrgb = convert_to_RGB(IM, cm, cmLIM)
% convert_to_RGB - converts any image to truecolor RGB using a specified colormap
% USAGE: IMrgb = convert_to_RGB(IM, cm, cmLIM)
% INPUTS:
%    IM    = the image [m x n]
%    cm    = the colormap [p x 3], optional; default = jet(256)
%    cmLIM = the data limits [min max] to be used in the color-mapping
%            optional; default = [min(IM) max(IM)]
% OUTPUTS:
%    IMrgb = the truecolor RGB image [m x n x 3]
% Based on ind2rgb from the Image Processing Toolbox
% EA Allen August 30, 2011
% eallen@mrn.org
% modified by wolf zinke, Sep 2014

%--------------------------------------------------------------------------
    if(nargin < 2)
        cm = jet(256);
    end
    nIND = size(cm,1);

    % set values within defined value range
    if(nargin < 3)
        cmLIM = [min(IM(:)) max(IM(:))];
    elseif(length(cmLIM) == 1)
        cmLIM = [-cmLIM cmLIM];
    end

    % normalize image data to fit the color map
    IM = IM - cmLIM(1);
    IM = IM / (cmLIM(2) - cmLIM(1));
    IM(IM<0) = 0;
    IM(IM>1) = 1;
    IM = round(IM*(nIND-1));

    IM = double(IM)+1;

    % define color maps
    r = zeros(size(IM));
    g = zeros(size(IM));
    b = zeros(size(IM));
    r(:) = cm(IM,1);
    g(:) = cm(IM,2);
    b(:) = cm(IM,3);

    % Fill in the r, g, and b channels
    IMrgb = zeros([size(IM),3]);
    IMrgb(:,:,1) = r;
    IMrgb(:,:,2) = g;
    IMrgb(:,:,3) = b;

%--------------------------------------------------------------------------
%% define color map
% Matlabs default jet colormap is a pain and not well suited for signed data.
function cmap = twowaycol(m)
% check https://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
% http://colorbrewer2.org/

    if ~nargin
        m = 64;
    end

    % GMT no green (NICE)
    colmat = [32,96,255; 32,159,255; 32,191,255; 0,207,255; 42,255,255; 85,255,255; ...
              127,255,255; 170,255,255; 255,255,84; 255,240,0; 255,191,0; 255,168,0; ...
              255,138,0; 255,112,0; 255,77,0; 255,0,0] ./ 255;
% 
%       sunsetred = [174, 208, 210, 237, 245, 249, 255, 255, 230, 180, 153, 119, 58, 0, 61];
%       sunsetgreen = [28, 50, 77, 135, 162, 189, 227, 250, 245, 221, 199, 183, 137, 139, 82];
%       sunsetblue = [62, 50, 62, 94, 117, 126, 170, 210, 254, 247, 236, 229, 201, 206, 161];
%       colmat = flipud([sunsetred(:),  sunsetgreen(:),  sunsetblue(:)]) ./ 255;

    % % % RdYlBu10 (GOOD)
    % colmat = flipud([0.6471      0 0.1490; 0.8431 0.1882 0.1529; 0.9569 0.4275 0.2627; ...
    %                  0.9922 0.6824 0.3804; 0.9961 0.8784 0.5647; 0.8784 0.9529 0.9725; ...
    %                  0.6706 0.8510 0.9137; 0.4549 0.6784 0.8196; 0.2706 0.4588 0.7059; ...
    %                  0.1922 0.2118 0.5843]);

    % % GMT_panoply (NICE)
%     colmat = [0.015686 0.054902 0.847059;0.125490 0.313725 1.000000;0.254902 0.588235 1.000000;0.427451 0.756863 1.000000; ...
%               0.525490 0.850980 1.000000;0.611765 0.933333 1.000000;0.686275 0.960784 1.000000;0.807843 1.000000 1.000000; ...
%               1.000000 0.996078 0.278431;1.000000 0.921569 0.000000;1.000000 0.768627 0.000000;1.000000 0.564706 0.000000; ...
%               1.000000 0.282353 0.000000;1.000000 0.000000 0.000000;0.835294 0.000000 0.000000;0.619608 0.000000 0.000000];
% 
    % % cmp_flux  (CLEAR)
%     colmat = [0 253 253;  8 222 253;  16 189 253;  24 157 253;  32 125 253;  40 93 253;  48 60 253;  85 85 253;  133 133 253;  181 181 253; ...
%     230 230 253;  253 230 230;  253 181 181;  253 133 133;  253 85 85;  253 60 48;  253 93 40;  253 125 32;  253 157 24;  253 189 16; ...
%     253 224 8;  253 253 0] ./ 255;

%     % lbmap (DECENT)
%     colmat = flipud([175  53  71; 216  82  88; 239 133 122; 245 177 139; 249 216 168; 242 238 197;
%                      216 236 241; 154 217 238; 68 199 239; 0 170 226; 0 116 188]/255);

    cmap = interp1(linspace(0,1,size(colmat,1)),colmat,linspace(0,1,m),'cubic');


%--------------------------------------------------------------------------
%% helper function to remove the white space added by Matlab to the image
% No clue how better do a work around to get rid of stupid Matlab behaviour
% thanks to Andy Bliss for putting this on File Excahnge
% http://www.mathworks.com/matlabcentral/fileexchange/20427-crop-whitespace-from-an-image
function crop(filename,append,margin)
%CROP gets rid of whitespace around an image
%   CROP(FILENAME,APPEND,MARGIN) is the full calling form. APPEND and
%      MARGIN are optional inputs.
%
%   CROP('filename.ext') crops the image in the file and saves it using
%   the original filename, overwriting the old image. The extension (ext)
%   can be anything IMREAD supports.
%
%   CROP(directory) crops all images in a directory.
%
%   If APPEND is 1, CROP saves the cropped image as 'filename_crop.ext'
%   in the same directory as the original.
%
%   MARGIN sets the margin width in pixels (default is 10).
%
%   Changes since version 1:
%       1. Now accepts directories as input.
%       2. Improved input checks and comments.
%       3. Handles transparency in .png files.
%       4. Just prior to saving version 1, I added the APPEND option.
%
%   Example:
%       crop('C:\MATLAB7\toolbox\matlab\demos\html\cruller_01.png',1)
%
%   See also: IMREAD, IMWRITE

%   Requirements: your FIND must allow the 'last' option (version 7+?)
%   Written by Andy Bliss Sept 8th, 2006. Revised May 31, 2012.

    %set the default margin width in pixels
    if nargin<3 || isempty(margin)
        margin=10; %15;
    end
    %default is to save with original filename
    if nargin<2 || isempty(append)
        append=0;
    end

    %Get image names
    if isstruct(filename) %assuming it is a struct as produced from DIR
        files=filename;
    elseif isdir(filename) %if the input is a directory, get all the image files from the directory
        currentdir=pwd;
        cd (filename)
        files=[dir('*.png'); dir('*.gif'); dir('*.bmp'); dir('*.jpg')];
    else %if it is a single file:
        files.name=filename;
    end

    %loop over all the files
    for n=1:length(files)
        filename=files(n).name;

        %get file info
        info=imfinfo(filename);

        %get the image
        if strcmp(info.Format,'png') %if it has transparent pixels
            T=imread(filename,'backgroundcolor',[1 1 1]); %backgroundcolor makes transparent pixels white, so they don't affect cropping.
        else
            T=imread(filename);
        end

        %sum the RGB values of the image
        xsum=sum(sum(T,3));
        ysum=sum(sum(T,3),2);
        % figure,plot(xsum),title('xsum'),xlabel('distance from left (pixels)'),ylabel('image intensity (big numbers are white)')
        % figure,plot(ysum),title('ysum'),xlabel('distance from top (pixels)'),ylabel('image intensity (big numbers are white)')

        %xsum will be equal to max(xsum) wherever there is a blank column in
        %   the image (rgb white is [255,255,255]). The left edge for the
        %   cropped image is found by looking for the first column in which
        %   xsum is less than max(xsum) and then subtracting the margin.
        %   Similar code for other edges.
        xleftedge=find(xsum<max(xsum),1,'first')-margin;
        if xleftedge<1
            xleftedge=1;
        end
        xrightedge=find(xsum<max(xsum),1,'last')+margin;
        if xrightedge>length(xsum)
            xrightedge=length(xsum);
        end
        ytopedge=find(ysum<max(ysum),1,'first')-margin;
        if ytopedge<1
            ytopedge=1;
        end
        ybottomedge=find(ysum<max(ysum),1,'last')+margin;
        if ybottomedge>length(ysum)
            ybottomedge=length(ysum);
        end

        %resave the image
        if append
            filename=[filename(1:end-4) '_crop' filename(end-3:end)];
        end
        imwrite(T(ytopedge:ybottomedge,xleftedge:xrightedge,:),filename)
    end
    %change back to calling directory, if necessary
    if exist('currentdir','var')
        cd(currentdir)
    end