function [out_img] = gaussian_filter(input_img, sigma)
    %% Filter setup
    % radius of affected pixels
    N = 5;
    ind = -floor(N/2) : floor(N/2);
    [X,Y] = meshgrid(ind,ind);
    % gauss function
    h = exp(-(X.^2 + Y.^2) / (2*sigma+sigma));
    h = h/ sum(h(:));

    % convert to column vector
    h = h(:);

    
    %% filter image
    I = im2double(input_img);
    I_pad = padarray(I, [floor(N/2) floor(N/2)]);
    C = im2col(I_pad, [N N], 'sliding');
    C_filter = sum(bsxfun(@times, C, h), 1);
    out_img = col2im(C_filter, [N N], size(I_pad), 'sliding');
end