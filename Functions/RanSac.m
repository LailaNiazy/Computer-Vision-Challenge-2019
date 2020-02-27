function [Korrespondenzen_robust]=RanSac(Korrespondenzen,num,K)

%% Vorbereitung
    %Parameter
    epsilon = 0.65;%0.5
    p=0.95;%0.5
    tolerance=5;%0.04;
    %num=8; % vorher k
    s=log(1-p)/log(1-(1-epsilon)^num); % vorher ^k
    largest_set_size=0;
    largest_set_dist=Inf;


%% Ausf�hren von RanSac
    i = 1;
    %while i <= s
    while largest_set_size < num %macht den Code sehr langsam
        i = i + 1;
        %choose k random columns
        x = randi(size(Korrespondenzen,2),1,num);
        columns = Korrespondenzen(:,x);
        [x1_pixel,x2_pixel,F] = achtpunktalgorithmus(columns,K);
        
        %calculate Sampson-Distance
        sd =  sampson_dist(F,x1_pixel, x2_pixel);
        a = (sd < tolerance).*sd;
        [row,col,v] = find(a);
        indices = sub2ind(size(a),row,col);
        set_size = length(v);
        dist = sum(v);
        if set_size > largest_set_size
           largest_set_size = set_size;
           Korrespondenzen_robust = Korrespondenzen(:,indices);
        elseif set_size == largest_set_size && dist < largest_set_dist
           largest_set_dist = dist;
           Korrespondenzen_robust = Korrespondenzen(:,indices);
        end
    end
%     for i = 1:size(Korrespondenzen,2)
%         % Schaetzung der essentiellen bzw fundamentalen Matrix und
%         % �bergeben von homogenen Koordinaten
%         neuMatrix = randperm(length(Korrespondenzen),k);
%         [x1_pixel,x2_pixel,EF] = achtpunktalgorithmus(Korrespondenzen(:,neuMatrix));%%%%%%%%%
% 
%         Consensus_Set = [];
%         dist = 0;
% 
%         % berechnung der sampson distance
%         sd = sampson_dist(EF, x1_pixel, x2_pixel);
% 
% 
%         % pruefung ob f�r Consensus_Set geeignet
%         for j = 1:size(sd,2)
%             if sd(j) < tolerance
% 
%                 Consensus_Set = [Consensus_Set, [x1_pixel(1:2,j); x2_pixel(1:2,j)]];
%                 dist = dist + sd;
%             end
%         end   
%         set_size = size(Consensus_Set,2);
%         
%         % pruefung ob neues set an Korrespondzen bedingungen f�r largest_set erfuellt
%         if set_size > largest_set_size
%             largest_set_size = set_size;
%             largest_set_dist = dist; 
%             Korrespondenzen_robust = Consensus_Set;
%         elseif set_size == largest_set_size
%                 if dist < largest_set_dist
%                     largest_set_size = set_size;
%                     largest_set_dist = dist; 
%                     Korrespondenzen_robust = Consensus_Set;
%                 end
%         end
% 
% 
%     end
    Korrespondenzen_robust;
end

%% Sampson Distanz
function d=sampson_dist(F, x1_pixel, x2_pixel)
    e=[0;   0;   1 ];
    e = dach(e);
    d=sum(x2_pixel.*(F*x1_pixel)).^2 ./(sum((e*F*x1_pixel).^2)+ sum((e*F'*x2_pixel).^2));
end
