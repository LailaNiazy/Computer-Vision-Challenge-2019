function [im3,disp3] = initialSynthesize(im0,im1,disp0,disp1,p)
%% Synthesize initial view and depth map
% initial synthesis with two rectified views and their corresponding disparity map
%% code
disp('----------------------------------starting with view synthesis...');
time0 = cputime;

%% 1st Step Initialization
vld0 = disp0~=0;

vld1 = disp1~=0;

[imh,imw,~] = size(im0);

dist0 = 1./disp0;
dist1 = 1./disp1;

dist_trd = 0.0001; %threshold value

sigfactor = 5;

im2 = zeros(imh,imw,3);
dist2 = 99 * ones(imh,imw);
im2w = zeros(imh,imw) + 0.00000001;
im3 = zeros(imh,imw,3);
dist3 = 99 * ones(imh,imw);
im3w = zeros(imh,imw) + 0.00000001;


%% 2nd Step: synthesis each view using the method 3D Warping
%  since the reference view im0 and im1 are rectified meaning they are located on the baseline: 
%  3D Warping is performed by shifting the pixels in the reference views horizontally by its scaled disparity 
%  
for i = 1:imh
    for j = 1:imw
        if (vld0(i,j) ~= 0)
            ncid = (j - p * disp0(i,j)); %index for the virtual view depending on p
            ncidl = floor(ncid);
            ncidr = ceil(ncid);
            wl = exp(-sigfactor * (ncid - ncidl)); 
            wr = exp(-sigfactor * (ncidr - ncid));
            if (ncidl > 0)
                if (abs(dist2(i,ncidl) - dist0(i,j)) < dist_trd)
                    im2(i,ncidl,:) = im2(i,ncidl,:) + im0(i,j,:) * wl;
                    im2w(i,ncidl,:) = im2w(i,ncidl,:)+ wl;
                elseif (dist2(i,ncidl) > dist0(i,j))
                    dist2(i,ncidl) = dist0(i,j);
                    im2(i,ncidl,:) = im0(i,j,:) * wl;
                    im2w(i,ncidl,:) = wl;
                end
            end
            if (ncidr > 0)
                if (abs(dist2(i,ncidr) - dist0(i,j)) < dist_trd)
                    im2(i,ncidr,:) = im2(i,ncidr,:) + im0(i,j,:) * wr;
                    im2w(i,ncidr,:) = im2w(i,ncidr,:)+ wr;
                elseif (dist2(i,ncidr) > dist0(i,j))
                    dist2(i,ncidr) = dist0(i,j);
                    im2(i,ncidr,:) = im0(i,j,:) * wr;
                    im2w(i,ncidr,:) = wr;
                end
            end
        end
    end
end
%%% warping done for the second reference image
for i = 1:imh
    for j = 1:imw
        if (vld1(i,j) ~= 0)
            ncid = (j + (1 - p) * disp1(i,j));
            ncidl = floor(ncid);
            ncidr = ceil(ncid);
            wl = exp(-sigfactor * (ncid - ncidl));
            wr = exp(-sigfactor * (ncidr - ncid));
            if (ncidl <= imw)
                if (abs(dist3(i,ncidl) - dist1(i,j)) < dist_trd)
                    im3(i,ncidl,:) = im3(i,ncidl,:) + im1(i,j,:) * wl;
                    im3w(i,ncidl,:) = im3w(i,ncidl,:)+ wl;
                elseif (dist3(i,ncidl) > dist1(i,j))
                    dist3(i,ncidl) = dist1(i,j);
                    im3(i,ncidl,:) = im1(i,j,:) * wl;
                    im3w(i,ncidl,:) = wl;
                end
            end
            if (ncidr <= imw)
                if (abs(dist3(i,ncidr) - dist1(i,j)) < dist_trd)
                    im3(i,ncidr,:) = im3(i,ncidr,:) + im1(i,j,:) * wr;
                    im3w(i,ncidr,:) = im3w(i,ncidr,:)+ wr;
                elseif (dist3(i,ncidr) > dist1(i,j))
                    dist3(i,ncidr) = dist1(i,j);
                    im3(i,ncidr,:) = im1(i,j,:) * wr;
                    im3w(i,ncidr,:) = wr;
                end
            end
        end
    end
end
nim2 = im2./repmat(im2w,1,1,3);
nim3 = im3./repmat(im3w,1,1,3);

%% 3rd Merge two views 
im3 = zeros(imh,imw,3); %virtual view
disp3 = min(dist2,dist3); %virtual disparity map

% only if both DepthMaps are under a specified threshold, both are blended
% using p, where the the closer reference image is assigned a higher weight
% otherwise it takes the pixel value of the reference image with the
% highest depth
for i = 1:imh
    for j = 1:imw
        if (abs(dist2(i,j) - dist3(i,j)) < dist_trd) 
            im3(i,j,:) = nim2(i,j,:)*(1-p) + nim3(i,j,:)*(p);
        elseif (dist2(i,j) > dist3(i,j))
            im3(i,j,:) = nim3(i,j,:);
        else
            im3(i,j,:) = nim2(i,j,:);
        end
    end
end

%% for displaying purposes
disp3 = disp3 + (-100) * (disp3 == 99);

time1 = cputime;
fprintf('----------------------------------initiate synthesize completed: %.2fs ... \n',time1 - time0);

end
