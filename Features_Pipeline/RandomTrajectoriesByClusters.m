function [] = RandomTrajectoriesByClusters(clusters,xys,...
    CellID, CellType)
% Input: clusters - numeric array with cluster of each cell
%        xys - cell array with each cell`s trajectory
%        CellID - numeric array with cells ids
%        CellType - categorical array with cells types
% Output: figures with 16 random cell trajectories from each cluster

% number of trajectories in each figure
Nc = 16;

% creating colormap for each trajectory
cidk = jet(length(clusters)-1);

%running over all clusters
for i=1:length(unique(clusters))
    % getting current clusters
    temp_clust = i;
    
    % find the cells consists the relevant cluster
    clust_idx = find(clusters == temp_clust);
    
    %new figure
    figure();
   
    if length(clust_idx) >= Nc
        random_idx = randperm(length(clust_idx),Nc);
        for k=1:1:Nc
            % getting sepcific cell trajectory
            xy=xys{clust_idx(random_idx(k))}; 
            
            % zero cell position at first frame.
            xy=xy-ones(size(xy(:,1)))*xy(1,:);  
            
            %plotting
            subplot(sqrt(round(Nc)),sqrt(round(Nc)),k)
            plot(xy(:,1),xy(:,2),'-','color', cidk(i,:), 'linewidth',1.5);
            title(strcat('Cluster:', {' '}, num2str(i), {' '}, 'CellID:',...
                {' '},num2str(CellID(clust_idx(random_idx(k)))), {' '},...
                'Cell Type:', {' '}, char(CellType(clust_idx(random_idx(k))))));
            set(gca,'xtick',[],'ytick',[]);
        end
    elseif length(clust_idx) < Nc
        for k=1:1:length(clust_idx)
            % getting sepcific cell trajectory
            xy=xys{clust_idx(k)}; 
            
            % zero cell position at first frame.  
            xy=xy-ones(size(xy(:,1)))*xy(1,:); 
            
            %plotting
            subplot(sqrt(round(Nc)),sqrt(round(Nc)),k)
            plot(xy(:,1),xy(:,2),'-','color', cidk(i,:), 'linewidth',1.5);
            title(strcat('Cluster', {' '}, num2str(i), {' '}, 'CellID',...
                {' '},num2str(CellID(clust_idx(k))), {' '}, 'Cell Type:', {' '},...
                char(CellType(clust_idx(k)))));
            set(gca,'xtick',[],'ytick',[]);
        end
    end
end
end

