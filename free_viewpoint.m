function output_image  = free_viewpoint(imageL, imageR, p,Ki, depthMap)
% This function generates an image from a virtual viewpoint between two
% real images. The output image has the same size as the input images.

%% Calculate correspondencies, essential matrix and R,T,lambda
display('start calculating correspondencies ...');
Korrespondenzen_robust = findKorrespondenzen(imageL,imageR,Ki);

display('finished calculating correspondencies ...');
display('start calculating rotation and translation ...');
[~,~,EF] = achtpunktalgorithmus(Korrespondenzen_robust,Ki);
[T,R,~,~] = euklidischeTransformation(EF,Ki,Korrespondenzen_robust);

display('finished calculating rotation and translation ...');
%% first import depthmap and rectify images

display('start rectifying ...');
[im0,im1,TL] = rectifyImageE(imageL,imageR,R,T,Ki);
display('loading DepthMap ...');
if depthMap == 1
    clear DepthMapL DepthMapR disparityFl disparityFr;
    load('DisparityMap/New_Depth_Map1.mat');
elseif depthMap == 2
    clear DepthMapL DepthMapR disparityFl disparityFr;
    load('DisparityMap/New_Depth_Map2.mat');
else
    %error('Could not load a valid depth map. Please make sure you added everything to the current path.');
    clear DepthMapL DepthMapR disparityFl disparityFr;
    fprintf('Start calculating disparitymap ... \n')
    [disparityFl, ~, ~, ~] = stereomatch(rgb_to_gray(im0), rgb_to_gray(im1), 7, 100);
    [disparityFr, ~, ~, ~] = stereomatch(rgb_to_gray(im1), rgb_to_gray(im0), 7 , 100);
end

display('finished loading ...');
display(size(disparityFl));




%% synthesizing new image
% remove ghost contour artifacts
[im0,im1] = removeGhostContour(im2double(im0),im2double(im1),disparityFl,disparityFr);
% synthesize depth map for new view and then the image  
[smooth_view,~] = initialSynthesize(im0,im1,disparityFl,disparityFr,p);

%% generating output_image

output_image = smooth_view;

end
