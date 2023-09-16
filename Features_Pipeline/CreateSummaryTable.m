%% Load variables
uniq_patients = unique(data.Patient);
uniq_clusters = unique(clusters);
uniq_celltype = unique(data.CellType);

%% Create Table
patient_cells_per_cluster = [];
cells_type_per_cluster = [];
cells_gender_per_cluster = [];

%patient cells per clusters
for i=1:length(uniq_clusters)
    temp_summary = summary(table(data.Patient(clusters == i)));
    patient_cells_per_cluster = [patient_cells_per_cluster; temp_summary.Var1.Counts'];
end

table1 = array2table(patient_cells_per_cluster,'VariableNames', temp_summary.Var1.Categories',...
    'RowNames', string(uniq_clusters));

%cells type per cluster
for i=1:length(uniq_clusters)
    temp_summary = summary(table(data.CellType(clusters == i)));
    cells_type_per_cluster = [cells_type_per_cluster; temp_summary.Var1.Counts'];
end

table2 = array2table(cells_type_per_cluster,'VariableNames', temp_summary.Var1.Categories',...
    'RowNames', string(uniq_clusters));

%cells gender per cluster
for i=1:length(uniq_clusters)
    temp_summary = summary(table(data.Gender(clusters == i)));
    cells_gender_per_cluster = [cells_gender_per_cluster; temp_summary.Var1.Counts'];
end

table3 = array2table(cells_gender_per_cluster,'VariableNames', temp_summary.Var1.Categories',...
    'RowNames', string(uniq_clusters));

Tright = join(table2,table3, 'Keys', 'Row');
ClustersTable = join(Tright,table1, 'Keys', 'Row');

