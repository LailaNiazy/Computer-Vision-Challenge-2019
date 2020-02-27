function Korrespondenzen_robust = findKorrespondenzen(imageL,imageR,Ki)
sigma=0.3;

grayL=rgb_to_gray(imageL);
grayR=rgb_to_gray(imageR);

grayL=gaussian_filter(grayL, sigma);
grayR=gaussian_filter(grayR, sigma);

%% harris parameters
% empirically determined parameters which find nice correspondencies
% as usual, these values would never have been found with logical thinking
% since they do not really make sense
min_dist=30;
tau=10e-4;
k=0.02;
segment_length=3;
N=400;
tile_size=[20 30];
do_plot = false;

HarrisL = harris_detektor(grayL, 'segment_length', segment_length, ...
    'k',k,'N',N,'tau',tau,'min_dist',min_dist,'tile_size',tile_size, ...
    'do_plot',do_plot);
HarrisR = harris_detektor(grayR, 'segment_length', segment_length, ...
    'k',k,'N',N,'tau',tau,'min_dist',min_dist,'tile_size',tile_size, ...
    'do_plot',do_plot);

%% parameters to match correspondencies
window_length=13;
minimum_correlation=0.9;

[Korrespondenzen] = korrespondenzen(grayL, grayR, HarrisL,HarrisR, ...
    window_length, minimum_correlation);

%% filter correspondencies to find relevant ones
Korrespondenzen_robust=[];

num = 8;
% try until a usable number of correspondencies is returned
while (size(Korrespondenzen_robust,2) < num)
    [Korrespondenzen_robust]=RanSac(Korrespondenzen,num,Ki);
end

if do_plot
    title = strcat( ...
        'Harriskorrespondenzen: min dist=',num2str(min_dist), ...
        '/segment length=',num2str(segment_length), ...
        '/tau=',num2str(tau),'/tile size=',num2str(tile_size), ...
        '/k=',num2str(k),'/N=',num2str(N));
    plot_corres(grayL, grayR, Korrespondenzen_robust,titel);
end

end
