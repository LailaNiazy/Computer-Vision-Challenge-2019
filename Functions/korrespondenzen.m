function [Korrespondenzen] = korrespondenzen(Im1, Im2, Mpt1, Mpt2,window_length,min_corr)
% finding and giving back correspondencies from the images

%% image preprocessing
tempMpt1 = [];
tempMpt2 = [];
cutoff = floor(window_length/2);
[r1,c1] = size(Im1);
cend1 = c1-cutoff;
rend1 = r1-cutoff;

%cutting features
for i = 1:size(Mpt1,2)
    if(Mpt1(1,i)>cutoff && (Mpt1(1,i)<=cend1) && Mpt1(2,i)>cutoff && (Mpt1(2,i)<=rend1))
        tempMpt1=[tempMpt1,Mpt1(:,i)];
    end
end

[r2,c2] = size(Im2);
cend2 = c2-cutoff;
rend2 = r2-cutoff;

for i = 1:size(Mpt2,2)
    if(Mpt2(1,i)>cutoff && (Mpt2(1,i)<=cend2) && Mpt2(2,i)>cutoff && (Mpt2(2,i)<=rend2))
        tempMpt2=[tempMpt2,Mpt2(:,i)];
    end
end

Mpt1=tempMpt1;
Mpt2=tempMpt2;

%% Normierung Fenster
Mat_feat_1 = norm_features(Im1,Mpt1,window_length);
Mat_feat_2 = norm_features(Im2,Mpt2,window_length);

%% NCC Brechnung
    % Anzahl Pixel
    N = window_length*window_length;
    NCC_matrix=double(zeros(size(Mat_feat_2,2),size(Mat_feat_1,2)));
    % Berechnung der Eintraege der NCC_matrix
    for i = 1:size(Mat_feat_1,2)
        for j = 1:size(Mat_feat_2,2)    
            NCC_matrix(j,i) = double((1/(N-1))*dot(Mat_feat_1(:,i)',Mat_feat_2(:,j)));
        end
    end
    % Eliminierung der Werte kleiner des Schwellwertes min_corr
    NCC_matrix = NCC_matrix.*(NCC_matrix>=min_corr);
    % Sortierung in absteigender Reihenfolge
    %NCC_matrix2 = NCC_matrix;
    [~,d] = sort(NCC_matrix(:), 'descend');
    % Speicherung in sorted_index
    nonzero = numel(find(NCC_matrix));
    sorted_index = d(1:nonzero,1);

%% Korrespondezen

% initialize the correspondence matrix and the size of NCC_matrix
    Korrespondenzen = [];
    [sz_row, sz_col] = size(NCC_matrix);
    for i = 1:numel(sorted_index)
            % bestimmung der Punkte y in Bild 1 und x in Bild 2
            [y, x] = ind2sub([sz_row, sz_col], sorted_index(i));
            % check if the value is 0, if no add the correspondence and
            % delete the column
            if NCC_matrix(y,x) == 0
            else
                Korrespondenzen = [Korrespondenzen,[Mpt1(:,x);Mpt2(:,y)]];
                NCC_matrix(:,x) = 0;
            end
    end
end

function Mat_feat=norm_features(I,Mpt,window_length)
    half_wl = floor(window_length/2);
    % pixel number in window
    n = window_length^2;
    
    % initialization
    Mat_feat = zeros(n,size(Mpt,2));
    
    for i = 1:size(Mpt,2)
        % get the window
        x = Mpt(1,i);
        y = Mpt(2,i);
        
        sy = y - half_wl;
        ey = y + half_wl;
        sx = x - half_wl;
        ex = x + half_wl;
        
        temp = double(I(sy:ey,sx:ex));
        
        % calculate average of temp
        average = sum(sum(temp))/n;

        % divide by standard deviation
        st_dev = std(temp(:));
        
        temp = (temp-average)/st_dev;
        
        Mat_feat(:,i) = temp(:);
    end
end