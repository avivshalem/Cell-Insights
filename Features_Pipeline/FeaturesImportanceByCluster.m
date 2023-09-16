function [out] = FeaturesImportanceByCluster(data, feature_names, clusters)
% Input: data - features numeric array where each column is a feature and
% each row relates to different cell
%       feature_names - name of each column feature in data
%        clusters - clusters of the data
% Output: out - struct with feature importance for each cluster

%initialize output
out.Cluster = [];
out.feature_importance = [];

for i=1:length(unique(clusters))
    
    % running over each of the clusters
    temp_clust = i;

    % creating binary classification vecotr
    binary_classes = clusters == temp_clust;
    
    % create random forest for identifying the importance of the feature to
    % the cluster
    temp_model = fitcensemble(data,binary_classes, 'Method',...
    'bag','NumLearningCycles', 100);

    try
        % find the importance
        importance = oobPermutedPredictorImportance(temp_model);
        % Sort by Feature Importance
        [~,idx] = sort(importance,'descend');
        
        % Feature importance is descending from left to right
        feature_names_sorted_by_importance = feature_names(idx);
%         i, feature_names_sorted_by_importance
        % Insert feature importance in output "out" struct
        out.Cluster = [out.Cluster; temp_clust];
        out.feature_importance = [out.feature_importance; feature_names_sorted_by_importance];  
    catch
        warning('Error computing importance of features. Skipping importance for current cluster.');
    end  
end

end

