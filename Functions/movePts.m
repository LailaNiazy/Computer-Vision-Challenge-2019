function newPts=movePts(oldPts,R,T)
    sz=size(oldPts);
    newPts=zeros(sz);
    for i=1:max(sz)
        newPts(:,i)=R'*oldPts(:,i)-R'*T;
    end
end