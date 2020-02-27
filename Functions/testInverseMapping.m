%% Testing inverse mapping
tic
% %% Load images
% imageL = imread('img/L2.JPG');
% imageR = imread('img/R2.JPG');
imageL = imread('img/L2.JPG');
imageR = imread('img/R2.JPG');
%imageL = imread('testBilder/imageL.png');
%imageR = imread('testBilder/imageR.png');

load('Kalibrierungsmatrix.mat');
K=K2_opt;

p = 0.5;

%% get depth map
load('DepthMaps/Unrect_Depth_Map2.mat')
%depthL = imread('testBilder/depthL.png');
%depthR = imread('testBilder/depthR.png');
%depthL = double(depthL);
%depthR = double(depthR);

%to use
DMr2l=DepthMapr_unrect;
DMl2r=DepthMapl_unrect;
Dispr2l=disparityFr_unrect;
Displ2r=disparityFl_unrect;

%% Preprocess DepthMap
%load('transform_for_rect.mat')
%TODO save JL 
% DepthMap=JL;

% imageL=imresize(imageL,size(DepthMapl_unrect));
% imageR=imresize(imageR,size(DepthMapl_unrect));

%% preprocessing
grayL=rgb_to_gray(imageL);
grayR=rgb_to_gray(imageR);

%% TODO find interesting points 
%Harris
%RANSAC
%to get this
%Korrespondenzen_robust=[569 301 839 334; 2006 647 2191 632; 2569 1214 2433 1233; 1396 1004 1384 1009; 786 1262 664 1260; 2912 1635 2788 1699; 785 648 1045 663; 2300 1018 2427 1021; 627 727 721 744; 1726 439 1682 429]';
Korrespondenzen_robust=findKorrespondenzen(imageL,imageR,K);
plot_corres(imageL,imageR,Korrespondenzen_robust,'Korrespondenzen mit Harris');
%% estimating camera movement
[x1,x2,EF] = achtpunktalgorithmus(Korrespondenzen_robust,K);
[T,R,lambda,x1] = euklidischeTransformation(EF,K,Korrespondenzen_robust,false);

final = buildImage(imageL,imageR,Displ2r,Dispr2l,R,T,p,K);
toc
%% actual inverse mapping
function final = buildImage(imageL,imageR,Displ2r,Dispr2l,R,T,p,K)
sl=inverseMapping(imageL,Displ2r,[0;0;0],R,T,p,K);
sr=inverseMapping(imageR,Dispr2l,-R'*T,scaleR(R,-1),scaleT(T,-1),1-p,K);

% our merging process
%final=p*sl+(1-p)*sr;
%% Merging process
% other merging process
[imh,imw,imd] = size(imageL);
dist_trd = 0.0001;
dist2=1./Displ2r;
dist3=1./Dispr2l;

final= merging(imh,imw,dist2,dist3,dist_trd,p,sl,sr);

figure
imshow(final);
title('final image');
end

function mergedImage = merging(imh,imw,dist2,dist3,dist_trd,theta,nim2,nim3)
im3 = zeros(imh,imw,3);

disp3 = min(dist2,dist3);
for rid = 1:imh
    for cid = 1:imw
        if (abs(dist2(rid,cid) - dist3(rid,cid)) < dist_trd)
            im3(rid,cid,:) = nim2(rid,cid,:)*(1-theta) + nim3(rid,cid,:)*(theta);
        elseif (dist2(rid,cid) > dist3(rid,cid))
            im3(rid,cid,:) = nim3(rid,cid,:);
        else
            im3(rid,cid,:) = nim2(rid,cid,:);
        end
    end
end
mergedImage = im3;
end

function shifted=inverseMapping(image,DepthMap,imageOrigin,R,T,p,K)
% R and T for virtual view
Rv=scaleR(R,p);
Tv=scaleT(T,p);

% Origins
% c1 = [0 0 0]
% c2 = -R'*T;
c3=Rv'*imageOrigin-Rv'*Tv;

% warp depth map: forward mapping
% https://pure.tue.nl/ws/files/3198275/200910785.pdf
% p 96

% lambda_2*d_2 = K_2*R_2*K_1^-1 Z_1w d_1 - K_2*R_2*C_2
% assumption: K is camera calibration matrix
% assumption: all camera calibrations match K_1=K_2=K
% assumption: from C1 to C2 we rotate by R, so R=R_2

% assumption: d1 is homogenous depth map coordinate

sdm=size(DepthMap);
warped_depth_map=zeros([sdm 3]);

%% STEP 1

fprintf("Forward mapping of depth map to virtual camera\n");
msdm=sdm(1)*sdm(2);

origin = K*Rv*c3;
factor = K*Rv/K;

for y=1:sdm(1)
    progress(sdm(1),50,y);
    for x=1:sdm(2)
        warped_depth_map(y,x,:)=factor*DepthMap(y,x)*[x y 1]' - origin;
    end
end

%% check by plot
if false
    figure
    wdm_vec=reshape(warped_depth_map,[msdm,3]);
    scatter3(wdm_vec(:,1),wdm_vec(:,2),wdm_vec(:,3));
end

%% STEP 2

% se = offsetstrel('ball',5,5);

if true
%     % MATLAB
%     % 3D ball structure
%     rad=100;
%     [x,y,z]=ndgrid(-rad:rad);
%     se = strel(sqrt(x.^2 + y.^2 + z.^2) <= rad);
%     
%     %dilated warped
%     dwdm=zeros([sdm 3]);
%     %eroded dilated warped
%     edwdm=zeros([sdm 3]);
%     %eroded eroded dilated warped
%     eedwdm=zeros([sdm 3]);
%     
%     % 1. dilate
%     fprintf("Dilation of forward mapping \n");
%     dwdm = imdilate(warped_depth_map,se);
%     

%     % 2. erosion
%     fprintf("1. erosion of forward mapping \n");
%     edwdm = imerode(dwdm,se);
% 
%     % 3. erosion (again)
%     fprintf("2. erosion of forward mapping \n");
%     eedwdm = imerode(edwdm,se);

    %% TODO %% TODO %% TODO %% TODO %% TODO %% TODO %% TODO %% TODO %%
    eedwdm = warped_depth_map;
    %% TODO %% TODO %% TODO %% TODO %% TODO %% TODO %% TODO %% TODO %%
    
else
    % OURS (Doesn't work, is for 2D only)
%     dwdm = dilation(warped_depth_map);
%     edwdm = erosion(edwdm);
%    eedwdm = erosion (dwdm);
end

% %% check by plot
% if true
%     figure
%     mvec=reshape(dwdm,[msdm,3]);
%     scatter3(mvec(:,1),mvec(:,2),mvec(:,3));
% 
%     figure
%     mvec=reshape(edwdm,[msdm,3]);
%     scatter3(mvec(:,1),mvec(:,2),mvec(:,3));
%     
%     figure
%     mvec=reshape(eedwdm,[msdm,3]);
%     scatter3(mvec(:,1),mvec(:,2),mvec(:,3));
% end

%% STEP 3
fprintf("Calculation of corresponding world points \n");
%dilated warped
wps=zeros([sdm 3]);
factor=inv(Rv)*inv(K);

for y=1:sdm(1)
    progress(sdm(1),50,y);
    for x=1:sdm(2)
        d=[eedwdm(y,x,1);eedwdm(y,x,2);eedwdm(y,x,3)];
        wps(y,x,:) = c3 + factor*d;
    end
end

%% check by plot
if false
    figure
    mvec=reshape(wps,[msdm,3]);
    scatter3(mvec(:,1),mvec(:,2),mvec(:,3));
end

%% STEP 4.1
fprintf("Projection on left source texture image\n");
psti=zeros([sdm 3]);

%projection matrix
PM=[K zeros(3,1)]*[R -R*imageOrigin; zeros(1,3) 1];
for y=1:sdm(1)
    progress(sdm(1),50,y);
    for x=1:sdm(2)
        vec=[wps(y,x,1) wps(y,x,2) wps(y,x,3) 1]';
        psti(y,x,:)=PM*vec;
    end
end


% %% check by plot
% if true
%     figure
%     mvec=reshape(psti,[msdm,3]);
%     scatter3(mvec(:,1),mvec(:,2),mvec(:,3));
% end

% %% try to color it
% imageL=grayL(2:end,101:(end-101));
% imageR=grayR(2:end,101:(end-101));
% 
% if true
%     figure
%     hold on
%     %size((1-p)*imageL(:,:)+p*imageR(:,:))
%     %size(asdf(:,:,1))
%     imageL=rand(size(psti));
%     h=surf(psti(:,:,1),psti(:,:,2),psti(:,:,3),imageL,'EdgeColor','interp');
%     set(h,'FaceColor','interp','FaceAlpha',0.5,'EdgeAlpha',0.1,'EdgeColor','interp');
% end

%% norm to homogenous coordinates
fprintf("Norm to homogenous coordinates\n");
hc=zeros([sdm 2]);

for y=1:sdm(1)    
    progress(sdm(1),50,y);
    for x=1:sdm(2)
        vec=[psti(y,x,1) psti(y,x,2)]'/psti(y,x,3);
        hc(y,x,:)=vec;
    end
end

%%
fprintf("Fill with color from projected point\n");
colored=zeros([sdm 3]);
for y=1:sdm(1)
    progress(sdm(1),50,y);
    for x=1:sdm(2)
        % [y x] to assume matlab notation
        src=[hc(y,x,2) hc(y,x,1)];
        % check if we land in picture
        
        % start with black
        color=[0 0 0];
        % change if we have a corresponding image point
        if ((src(1)>=1) & (src(1)<=sdm(1)) & (src(2)>=1) & (src(2)<=sdm(2)))
            % grab color from four surrounding
            fy=floor(src(1));
            cy=ceil(src(1));
            fx=floor(src(2));
            cx=ceil(src(2));
            
            color=double((image(fy,fx,:)+image(fy,cx,:)+image(cy,fx,:)+image(cy,cx,:))/4);
        end
        colored(y,x,:)=color;
    end
end

% norm colors
colored=double(colored/max(max(max(colored(:,:,:)))));

%%
figure
imshow(colored);
title('Colored image');
shifted=colored;
end
