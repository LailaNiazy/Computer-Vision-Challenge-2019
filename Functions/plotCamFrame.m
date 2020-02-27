function plotCamFrame(cf,name,color)
    plot3([cf(1,:) cf(1,1)],[cf(2,:) cf(2,1)],[cf(3,:) cf(3,1)],color)
    text(cf(1,end),cf(2,end),cf(3,end),name,'Color',color);
end