function [T,R,lambda,x1] = euklidischeTransformation(EF,K,Korrespondenzen)
    [T1, R1, T2, R2, ~, ~]=TR_aus_E(EF);
    
    [T_cell, R_cell, d_cell, x1, x2] = tr_calib(T1, T2, R1, R2, Korrespondenzen, K);
    [T, R, lambda, ~, ~] = find_tr(T_cell, R_cell, d_cell, Korrespondenzen, x1, x2);
    lambda = abs(lambda);
    
    if false
        [P1, ~, ~] = showCamFramePosition(R, T, lambda, x1);
        [~, ~] = rueckprojektion(Korrespondenzen, P1, image, T, R, K);
    end
end

function [T1, R1, T2, R2, U, V]=TR_aus_E(E)
    % Diese Funktion berechnet die moeglichen Werte fuer T und R
    % aus der Essentiellen Matrix
    
    [U,S,V]=svd(E);
    
    if (det(U) ~= 1)
       % flip determinante by doing E=U*S*V'=U*E*S*V'
       % E=diag([1 1 -1])^2;
       % one half vanishes in S which is not changed, the other half modifies U
       U=U*diag([1 1 -1]);
    end
    
    if (det(V) ~= 1)
       % flip determinante by doing E=U*S*V'=U*S*E*V'
       % E=diag([1 1 -1])^2;
       % one half vanishes in S which is not changed, the other half modifies V
       % special care must be taken because of the transposing
       V=[diag([1 1 -1])*V']';
    end
        
    %Lecture 3.4 slide 4
    %R_z(+pi/2)
    R_z1=[0 -1 0; 1 0 0; 0 0 1];
    %R_z(-pi/2)
    R_z2=[0 1 0; -1 0 0; 0 0 1];

    R1 = U*R_z1'*V';
    R2 = U*R_z2'*V';
    
    T1dach = U*R_z1*S*U';
    T2dach = U*R_z2*S*U';
    
    T1=dedach(T1dach);
    T2=dedach(T2dach);
end

function w = dedach(W)
    % Diese Funktion implementiert die Umkehrung des ^-Operator.
    % Sie wandelt eine schiefsymmetrische Matrix in einen 3-Komponenten Vektor
    w=zeros(3,1);
    w(3)=W(2,1);
    w(2)=W(1,3);
    w(1)=W(3,2);
end

function [T_cell, R_cell, d_cell, x1, x2] = tr_calib(T1, T2, R1, R2, Korrespondenzen, K)
    %% provides several possibilities for R and T  
    % Preparation
    size_Korr = size(Korrespondenzen,2);
    
    x1=[Korrespondenzen(1:2,:);ones(1,size_Korr)];
    x2=[Korrespondenzen(3:4,:);ones(1,size_Korr)];
    
    % calibration of coordinates
    for i=1:size_Korr
        x1(:,i) = K\x1(:,i);
        x2(:,i) = K\x2(:,i);
    end
    
    % build possible combinations of T and R
    T_cell = {T1 T2 T1 T2};
    R_cell = {R1 R1 R2 R2};
    
    tmp = zeros(size_Korr,2);
    d_cell = {tmp tmp tmp tmp};
end

function [T, R, lambda, M1, M2] = find_tr(T_cell, R_cell, d_cell, Korrespondenzen, x1, x2)
    %% find best combination of T and R
    n=max(size(Korrespondenzen));
    M1=zeros(n*3,n+1);
    M2=M1;
    
    pos_count=0;
    best_i=0;

    for i=1:max(size(T_cell))
        M1=zeros(n*3,n+1);
        M2=M1;
        for j=1:n
            dx1=dach(x1(:,j)');
            dx2=dach(x2(:,j)');
            
            M1((j-1)*3+1:j*3,j)=dx2*R_cell{i}*x1(:,j);
            M1((j-1)*3+1:j*3,end)=dx2*T_cell{i};
            
            M2((j-1)*3+1:j*3,j)=dx1*R_cell{i}'*x2(:,j);
            M2((j-1)*3+1:j*3,end)=-dx1*R_cell{i}'*T_cell{i};
        end
        
        %Kapitel 5-1 - 5-3 Reconstruction from two calibrated views
        %page 125 (18 in pdf)
        
        %The linear least squares estimate of lambda is 
        %the eigenvector of M'M that corresponds to its smallest eigenvalue
        [~,~,V1]=svd(M1,0);
        [~,~,V2]=svd(M2,0);
        
        %The right-singular vectors of M (the colummns of V) are a set of orthonormal eigenvectors of M'∗M. (Wikipedia, svd)
        %Since the singular values are listed in descending order, lets assume the last column is the desired.
        lambda1=V1(:,end);
        lambda2=V2(:,end);
        
        d1=lambda1/lambda1(end);
        d2=lambda2/lambda2(end);
        
        d_cell{i}=[d1(1:end-1) d2(1:end-1)];
        
        curr_pos_count=sum(sum(d_cell{i}>0));
        
        if curr_pos_count > pos_count
            pos_count=curr_pos_count;
            best_i=i;
        end
    end
    
    T=T_cell{best_i};
    R=R_cell{best_i};
    lambda=d_cell{best_i};
end

function [P1, camC1, camC2] = showCamFramePosition(R,T,lambda, x1) 
    %% Visualizes the positions of the Cam frames and world points in world space
    % Annahme C1 sitzt im Ursprung. Koordinaten der Punkte sind auf dem Kameraframe
    % Also gibt x1 den Richtungsvektor des Punktes
    % lambda die dazugehörige Tiefeninfo
    
    P1=x1.*lambda(:,1)';
    
    camC1=[-0.2 0.2 0.2 -0.2;0.2 0.2 -0.2 -0.2;1 1 1 1];
    camC2=movePts(camC1,R,T);
    
    figure('Name','Rekonstruktion der Kameras')
    hold on;
    
    % plot reference points
    scatter3(P1(1,:),P1(2,:),P1(3,:),'k.');
    for i=1:max(size(P1))
        text(P1(1,i),P1(2,i),P1(3,i),num2str(i));
    end
    
    plotCamFrame(camC1,'Cam1','b');
    plotCamFrame(camC2,'Cam2','r');
    
    campos([43 -22 -87]);
    camup([0 -1 0]);
    
    style equal
    xlim auto; ylim auto; zlim auto;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on;
end

function [repro_error, x2_repro] = rueckprojektion(Korrespondenzen, P1, image, T, R, K)
    % Diese Funktion berechnet den mittleren Rueckprojektionsfehler der 
    % Weltkooridnaten P1 aus Bild 1 im Cameraframe 2 und stellt die 
    % korrekten Merkmalskoordinaten sowie die rueckprojezierten grafisch dar.
    
    % Punkte in P1 in Cameraframe 2 umrechnen
    P2=R*P1+T;
    
    % homogene koordinaten x2
    % P2 sind die Weltkoordinaten gegenüber Cam2 als Ursprung
    % wir haben die Tiefeninfos.
    % für x2 projezieren wir sie auf das Camframe, indem wir sie durch die Tiefeninfo teilen, und damit alle Punkte auf diese Ebene projezieren
    
    x2_repro = K*P2./P2(3,:);
    N=max(size(P2));
    repro_error=0;
    
    for i=1:N
        repro_error = repro_error + norm(Korrespondenzen(3:4,i)-x2_repro(1:2,i));
    end
    
    repro_error = repro_error/N;
    
    figure('Name','Rekonstruktion der Kamera 2')
    imshow(image);
    hold on;
    % draw both 
    plot(x2_repro(1,:),x2_repro(2,:),'r.');
    plot(Korrespondenzen(3,:),Korrespondenzen(4,:),'g.');

    for i=1:N
    % add numbers 
       txt=num2str(i);
       text(x2_repro(1,i),x2_repro(2,i),txt,'Color','r');
       text(Korrespondenzen(3,i),Korrespondenzen(4,i),txt,'Color','g');
    end
end
