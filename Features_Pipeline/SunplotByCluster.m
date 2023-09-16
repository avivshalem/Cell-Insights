function [] = SunplotByCluster(xys,clusters,maxtraj,saveflag, aviv, aviv2, loop, layer)
% Input: xys - Cell array consists of all cells trajectories.
%        clusters - clusters array. each row is the cluster of a different
%        cell.
%        maxtraj - maximum number of trajecotries in each sunplot
%        saveflag - 1 to save plot, 0 to not.
% Output: Sunplot by Clusters figure


% creating different color for each trajctory in each clusters

for i=1:length(unique(clusters))
    temp_clust = i;
    if aviv
        if aviv2
        name = unique(clusters);
        temp_clust = name(i);
        end
        clust_idx = find(clusters == temp_clust,maxtraj);
    else
        if aviv2
        name = unique(clusters);
        temp_clust = name(i);
        end
        clust_idx = find(clusters == temp_clust,maxtraj);
    end
    
    figure();
    cidk = jet(length(clust_idx)); 
    for k=1:1:length(clust_idx)
         xy=xys{clust_idx(k)};         
         xy=xy-ones(size(xy(:,1)))*xy(1,:); % zero cell position at first frame. 
         plot(xy(:,1),xy(:,2),'-','color',[cidk(k,:) 0.3],'linewidth',0.5); hold on; 
    end
    if aviv2
        title(strcat('Patient', {' '}, string(name(i))));
    else
        title(strcat('Cluster', {' '}, num2str(i)));
    end
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-350 350])
    xlim([-350 350])
%     
    if saveflag == 2
        if aviv2
            saveas(gcf,strcat('C:\Cell_Migration_MATLAB\Pipeline\Temp\', layer, '\', num2str(loop),'_SunPatTrans', num2str(i), '.png'));
        else
             saveas(gcf,strcat('C:\Cell_Migration_MATLAB\Pipeline\Temp\', layer, '\', num2str(loop),'_SunClustTrans', num2str(i), '.png'));
        end
    end
    if saveflag == 1
        mkdir('Graphs','Sunplots');
        addpath(strcat('Graphs/','Sunplots'));
        filename = char(strcat(pwd, '\Graphs\Sunplots', '\Sunplot of Cluster',{' '}, num2str(temp_clust),'.png'));
        saveas(gcf,filename);
    end
end
end

