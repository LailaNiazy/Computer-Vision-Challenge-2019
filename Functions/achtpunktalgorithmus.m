function[x1_pixel,x2_pixel,EF] = achtpunktalgorithmus(Korrespondenzen,K)
%% Preparation
% homogeneous Koordinaten
    x1_pixel=[Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen(1:2,:),2))];
    x2_pixel=[Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen(3:4,:),2))];

%% generate A and solve A*b=0
    
    % if possible, calibrate pixels
    if exist('K','var')
       x1_pixel=K\x1_pixel;
       x2_pixel=K\x2_pixel;
    end
    
    % do actual math
    A=zeros(size(x2_pixel,2),size(x1_pixel,1)*size(x2_pixel,1));
    for i=1:size(x2_pixel,2)
        A(i,:)=kron(x1_pixel(:,i),x2_pixel(:,i))';
    end
    
    % Solve
    [~,~,V]=svd(A);
 %% Projection onto next essential/fundamental matrix
 
    % solve min(E-G)
    G=reshape(V(:,9),[3,3]);

    [U,S,N]=svd(G);
    
    % fix singular values
    if exist('K','var')
        S = diag([1 1 0]);
        EF=U*S*N';
    else
        S(3,3)=0;
        EF=U*S*N';
    end
end
