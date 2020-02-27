%% plotting
function plot_corres(i1,i2,c,Titel)
        % This function overlays two images i1, i2 and plots and connects the given
        % correspondences c found in them 
        figure('Name',Titel);
        % plot Images
        imshow(i1);
        title(Titel);
        hold on;
        i2=imshow(i2);
        % make second transparent
        set(i2, 'AlphaData', 0.5);
        
        % plot correspondencies
        color = rand(size(c,2),3);
        
        scatter(c(1,:),c(2,:),3,color);
        scatter(c(3,:),c(4,:),3,color);
        % connect them
        for i=1:size(c,2)
            h = plot(c([1 3],i),c([2 4],i),'color','w');
            set(h(:),'linewidth',2);
        end
        hold off
end
