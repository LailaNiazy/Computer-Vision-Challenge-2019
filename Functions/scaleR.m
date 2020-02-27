function scaledR=scaleR(R,p)
    %drehachse finden
    a=R+R'-(trace(R)-1)*diag([1 1 1]);
    axis=a(:,1)/norm(a(:,1));
%     r=rand(3,1);
    % sichergehen, dass r nicht linear abhaengig ist

%     while (cross(axis,r)==zeros(3,1))
%         r=rand(3,1);
%     end
    r=(1/sqrt(axis(1)^2+axis(2)^2))*[-axis(2);axis(1);0];
    %vektor finden der normal auf a steht
    n=cross(axis,r);
 
    %vektor rotieren
    nrot = R*n;
    
%     axang = rotm2axang(R);
%     axang
%     phimat=rotm2eul(R);
%     phimat
%     newR=eul2rotm(phimat);
%     newR
%     R
    %winkel finden
    phi = acos(dot(n,nrot)/(norm(n)*norm(nrot)));
    
    newphi=p*phi;
    if(newphi<0)
        newphi=newphi+2*pi;
    end
    
    vers=1-cos(newphi)
    c=cos(newphi);
    s=sin(newphi);
	scaledR=[
        (axis(1)^2)*vers+c axis(1)*axis(2)*vers-axis(3)*s axis(1)*axis(3)*vers+axis(2)*s; 
        axis(1)*axis(2)*vers+axis(3)*s (axis(2)^2)*vers+c axis(2)*axis(3)*vers-axis(1)*s; 
        axis(1)*axis(3)*vers-axis(2)*s axis(2)*axis(3)*vers+axis(1)*s (axis(3)^2)*vers+c];
end
