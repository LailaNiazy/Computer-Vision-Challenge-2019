%% Computer Vision Challenge

% Groupnumber:
group_number = 58;

% Groupmembers:
members = {'Dominik Thoma', 'Falk Schoenfeld', 'Kathrin Khadra', 'Laila Niazy', 'Zuhra Amiri'};
 
% Email-Adress (from Moodle!):
mail = {'ga65qil@mytum.de', 'dominik.thoma@tum.de', 'falk.schoenfeld@tum.de', 'ga96tum@mytum.de', 'ga96mev@mytum.de','ga96cuh@mytum.de'};

%% add paths
addpath(genpath('RectifKitE'));
%addpath(genpath('RectifKitU'));
addpath(genpath('Functions'));
addpath(genpath('DisparityMap'));
addpath(genpath('Rendering_view_synthesis'));
%% Load images
imageL = imread('img/L2.JPG');
imageR = imread('img/R2.JPG');
p = 0.2;
load('Kalibrierungsmatrix.mat');

Ki = K2_opt;
DepthMap = 2;

  
%% Free Viewpoint -> normal execution
tic
output_image = free_viewpoint(imageL, imageR, p, Ki, DepthMap);
elapsed_time = toc;

figure()
imshow(output_image);
title('Virtual View 2 mit p = 0.2')
%3D Reconstruction from group 58');
display('Performance:');
display(elapsed_time);

%% Free Viewpoint -> execution with GUI
%guiTest();

%% Clear Workspace
%clear 