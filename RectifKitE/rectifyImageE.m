
function [JL,JR,TL] = rectifyImageE(IL,IR,R,T,K)
% Rectification with calibration data

% This function reads a (stereo) pair of images and respective camera matrices
% (PPMs) from files and rectify them. It outputs on files the two rectified
% images in PNG format. It reads  RGB images in PNG format.
%
% The bounding box and the transformation that has been applied
% are saved in the PNG metadata

%         Andrea Fusiello, 2007 (andrea.fusiello@univr.it)

% These points (in the left image) are used only for displaying the
% correspoding epipolar lines in the right image. They are not used in the
% rectification, therefore ml can be set to [0;0]


% read camera matrices pml and pmr
%load(['data/' img_base '_cam'])
% The user who wants to rectify his own images must provide 
% the two camra perspective projection matrices  pml, pmr here

%right origin
%left origin is 0 0 0 by definition
OL=zeros(3,1);
OR=R'*OL-R'*T;

%%
% At this point ml, pml and pmr are set.
pml=K*[R -T];
pmr=K*[R T];

% Epipolar geometry
[F,epil,epir] = fund(pml,pmr);
%ml=[2500 1900 723;533 300 768];

% --------------------  RECTIFICATION

disp('---------------------------------- rectifying...')

%  rectification without centeriing
[TL,TR,pml1,pmr1] = rectify(pml,pmr);

% centering LEFT image
p = [size(IL,1)/2; size(IL,2)/2; 1];
px = TL * p;
dL = p(1:2) - px(1:2)./px(3) ;

% centering RIGHT image
p = [size(IR,1)/2; size(IR,2)/2; 1];
px = TR * p;
dR = p(1:2) - px(1:2)./px(3) ;

% vertical diplacement must be the same
dL(2) = dR(2);

%  rectification with centering
[TL,TR,pml1,pmr1] = rectify(pml,pmr,dL,dR);

disp('---------------------------------- warping...')

% find the smallest bb containining both images
bb = mcbb(size(IL),size(IR), TL, TR);

% warp RGB channels,
for c = 1:3

    % this imwarp is NOT the matlab toolbox imwarp
    % Warp LEFT
    [JL(:,:,c),bbL,alphaL] = imwarp(IL(:,:,c), TL, 'bilinear', bb);

    % Warp RIGHT
    [JR(:,:,c),bbR,alphaR] = imwarp(IR(:,:,c), TR, 'bilinear', bb);

end

end

