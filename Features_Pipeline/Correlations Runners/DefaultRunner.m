%%% Correlations
%%%

load('\NatureAviv\HGPS_Data\DataFileWithHGPS.mat')
load('\NatureAviv\HGPS_Data\xys_full_with_HGPS.mat')
load('\NatureAviv\HGPS_Data\xys_TRJs_with_HGPS.mat')
load('\NatureAviv\Groups\GI_ONE_DIM_Merged_4_Use.mat')
load('\NatureAviv\Groups\GI_Good_7_Use.mat')
% toSever = (DataFile.CAP > 113) & (DataFile.CAP < 116);
% DataFile.CellType(toSever) = 'Severe';

DataFileWell = DataFileWithHGPS;
toSever = (DataFileWell.CAP > 113) & (DataFileWell.CAP < 500);
DataFileWell.CellType(toSever) = "'Severe'";
data = DataFileWell;

% find Wells
feature_vectors = [data.Dp, data.Dtot, data.Pp, data.Psi,...
    data.Pnp, data.Dnp, data.MSD10, ...
    data.Sp, data.Snp];

% TF = [];
% for feature=1:size(feature_vectors,2)
%     TF = [TF isoutlier(feature_vectors(:,feature))];
% end
% TF = isoutlier(feature_vectors, "median");
% TFRemoved = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);

WellsUniques = unique(data.CAPPatient2);
% CAPUniques = unique(data.CAP);
meanList = [];
meanStd = [];
for i=1:length(WellsUniques)
    currentData = find(data.CAPPatient2 == WellsUniques(i));
    currentFeatures = feature_vectors(currentData,7);

    meanList(1+end) = mean(currentFeatures);
    meanStd(1+end) = std(currentFeatures);
end
% scatter(WellsUniques, meanList)
wellListHC = WellsUniques(1:67);
wellListPremanifest = WellsUniques(138:187);
wellListMild = [WellsUniques(188:190); WellsUniques(68:84)];
wellListSevere = WellsUniques(85:137);
wellListHGPS = WellsUniques(191:end);
meanListHC = meanList(1:67);
meanListPremanifest = meanList(138:187);
meanListMild = [meanList(188:190) meanList(68:84)];
meanListSevere = meanList(85:137);
meanListHGPS = meanList(191:end);

wellCheckHC = wellListHC(meanListHC < 3300);
wellCheckPremanifest = wellListPremanifest(meanListPremanifest < 3200);
wellCheckMild = wellListMild(meanListMild < 0);
wellCheckSevere = wellListSevere(meanListSevere > 2100);
wellCheckHGPS = wellListHGPS(meanListHGPS > 1800);
%4
% wellCheckHC = wellListHC(meanListHC < 3700);
% wellCheckPremanifest = wellListPremanifest(meanListPremanifest < 3700);
% wellCheckMild = wellListMild(meanListMild > 3700 | meanListMild < 1700);
% wellCheckSevere = wellListSevere(meanListSevere > 1900);
% wellCheckHGPS = wellListHGPS(meanListHGPS > 1700);
%3
% wellCheckHC = wellListHC(meanListHC < 3700);
% wellCheckPremanifest = wellListPremanifest(meanListPremanifest < 3600);
% wellCheckMild = wellListMild(meanListMild < 0);
% wellCheckSevere = wellListSevere(meanListSevere > 1900);
% wellCheckHGPS = wellListHGPS(meanListHGPS > 1700);
%2
% wellCheckHC = wellListHC(meanListHC < 3700);
% wellCheckPremanifest = wellListPremanifest(meanListPremanifest < 3400);
% wellCheckMild = wellListMild(meanListMild < 0);
% wellCheckSevere = wellListSevere(meanListSevere > 1900);
% wellCheckHGPS = wellListHGPS(meanListHGPS > 1700);
% wellCheckHC = wellListHC(meanListHC < 3400);
%1
% wellCheckPremanifest = wellListPremanifest(meanListPremanifest < 3100);
% wellCheckMild = wellListMild(meanListMild < 0);
% wellCheckSevere = wellListSevere(meanListSevere > 2200);
% wellCheckHGPS = wellListHGPS(meanListHGPS > 1800);

disp(size(wellListHC,1) - size(wellCheckHC,1));
disp(size(wellListPremanifest,1) - size(wellCheckPremanifest,1));
disp(size(wellListMild,1) - size(wellCheckMild,1));
disp(size(wellListSevere,1) - size(wellCheckSevere,1));
disp(size(wellListHGPS,1) - size(wellCheckHGPS,1));
disp('************************');

%
% combinedIndex = [1:1:size(data.Dp,1)]';
% combinedIndex = and(~contains(string(data.Patient), '4709'), ~contains(string(data.Patient), '6274'));
% combinedIndex = ~contains(string(data.CAPPatient2), string(wellCheckHC)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckPremanifest)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckMild)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckSevere)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckHGPS)) & ...
%     ~contains(string(data.Patient), '4709') & ...
%     ~contains(string(data.Patient), '6274') & ...
%     TFRemoved;
% combinedIndex = ~contains(string(data.CAPPatient2), string(wellCheckHC)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckPremanifest)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckMild)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckSevere)) & ...
%     ~contains(string(data.CAPPatient2), string(wellCheckHGPS));
combinedIndex = ~contains(string(data.CAPPatient2), string(wellCheckHC)) & ...
    ~contains(string(data.CAPPatient2), string(wellCheckPremanifest)) & ...
    ~contains(string(data.CAPPatient2), string(wellCheckMild)) & ...
    ~contains(string(data.CAPPatient2), string(wellCheckSevere)) & ...
    ~contains(string(data.CAPPatient2), string(wellCheckHGPS)) & ...
    ~contains(string(data.Patient), '4709') & ...
    ~contains(string(data.Patient), '6274');
% combinedIndex = TFRemoved;

uniques = unique(data.CellType);
AllCountHC = find(data.CellType(combinedIndex) == uniques(1,:));
AllCountHGPS = find(data.CellType(combinedIndex) == uniques(2,:));
AllCountMild = find(data.CellType(combinedIndex) == uniques(3,:));
AllCountPremanifest = find(data.CellType(combinedIndex) == uniques(4,:));
AllCountSevere = find(data.CellType(combinedIndex) == uniques(5,:));
combinedIndexCached = combinedIndex;

AllCountHCAll = find(data.CellType == uniques(1,:));
AllCountHGPSAll = find(data.CellType == uniques(2,:));
AllCountMildAll = find(data.CellType == uniques(3,:));
AllCountPremanifestAll = find(data.CellType == uniques(4,:));
AllCountSevereAll = find(data.CellType == uniques(5,:));
%
close all
combinedIndex = combinedIndexCached;
% Correlations Loop Fixed X AXIS
folderCurrent = "Results\Regular";
mkdir(folderCurrent);
feature_vectors = [data.Dp, data.Dtot, data.Pp, data.Psi,...
    data.Pnp, data.Dnp, data.MSD10, ...
    data.Sp, data.Snp];
% feature_vectors = [data.Dp(combinedIndex), data.Dtot(combinedIndex), data.Pp(combinedIndex), data.Psi(combinedIndex),...
%     data.Pnp(combinedIndex), data.Dnp(combinedIndex), data.MSD10(combinedIndex), ...
%     data.Sp(combinedIndex), data.Snp(combinedIndex)];
% feature_names = [{'Dp'}, {'Dtot'},{'Pp'}, {'Psi'}, {'Pnp'}, {'Dnp'},...
%     {'MSD10'}, {'Sp'}, {'Snp'}];
TF = isoutlier(feature_vectors, "quartiles", ThresholdFactor=1.5);
TFRemoved = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
% num = 8.8; 1
% num = 6; 2
% num = 5.5; 3
% num = 2.75; 4 mean
num = 3.3777;

FactorHC = num;
FactorPre = num;
FactorMild = num;
FactorSevere = num;
FactorHGPS = num;

combinedIndexRemoved = false(size(combinedIndex));
currentCombinedIndex = false(size(combinedIndex));
currentCombinedIndex(AllCountHCAll) = combinedIndex(AllCountHCAll);
TF = isoutlier(feature_vectors(currentCombinedIndex,:), "mean", ThresholdFactor=FactorHC);
TFRemovedHC = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
TFIndex = 1;
for myIndex=1:length(currentCombinedIndex)
    if currentCombinedIndex(myIndex)
        if TFRemovedHC(TFIndex)
            if combinedIndexRemoved(myIndex)
                disp("******** ERROR ***********");
            end
            combinedIndexRemoved(myIndex) = true;
        end
        TFIndex = TFIndex + 1;
    end
end
currentCombinedIndex = false(size(combinedIndex));
currentCombinedIndex(AllCountPremanifestAll) = combinedIndex(AllCountPremanifestAll);
TF = isoutlier(feature_vectors(currentCombinedIndex,:), "mean", ThresholdFactor=FactorPre);
TFRemovedPremanifest = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
TFIndex = 1;
for myIndex=1:length(currentCombinedIndex)
    if currentCombinedIndex(myIndex)
        if TFRemovedPremanifest(TFIndex)
            if combinedIndexRemoved(myIndex)
                disp("******** ERROR ***********");
            end
            combinedIndexRemoved(myIndex) = true;
        end
        TFIndex = TFIndex + 1;
    end
end
currentCombinedIndex = false(size(combinedIndex));
currentCombinedIndex(AllCountMildAll) = combinedIndex(AllCountMildAll);
TF = isoutlier(feature_vectors(currentCombinedIndex,:), "mean", ThresholdFactor=FactorMild);
TFRemovedMild = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
TFIndex = 1;
for myIndex=1:length(currentCombinedIndex)
    if currentCombinedIndex(myIndex)
        if TFRemovedMild(TFIndex)
            if combinedIndexRemoved(myIndex)
                disp("******** ERROR ***********");
            end
            combinedIndexRemoved(myIndex) = true;
        end
        TFIndex = TFIndex + 1;
    end
end
currentCombinedIndex = false(size(combinedIndex));
currentCombinedIndex(AllCountSevereAll) = combinedIndex(AllCountSevereAll);
TF = isoutlier(feature_vectors(currentCombinedIndex,:), "mean", ThresholdFactor=FactorSevere);
TFRemovedSevere = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
TFIndex = 1;
for myIndex=1:length(currentCombinedIndex)
    if currentCombinedIndex(myIndex)
        if TFRemovedSevere(TFIndex)
            if combinedIndexRemoved(myIndex)
                disp("******** ERROR ***********");
            end
            combinedIndexRemoved(myIndex) = true;
        end
        TFIndex = TFIndex + 1;
    end
end
currentCombinedIndex = false(size(combinedIndex));
currentCombinedIndex(AllCountHGPSAll) = combinedIndex(AllCountHGPSAll);
TF = isoutlier(feature_vectors(currentCombinedIndex,:), "mean", ThresholdFactor=FactorHGPS);
TFRemovedHGPS = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
TFIndex = 1;
for myIndex=1:length(currentCombinedIndex)
    if currentCombinedIndex(myIndex)
        if TFRemovedHGPS(TFIndex)
            if combinedIndexRemoved(myIndex)
                disp("******** ERROR ***********");
            end
            combinedIndexRemoved(myIndex) = true;
        end
        TFIndex = TFIndex + 1;
    end
end
combinedIndex = logical(combinedIndexRemoved);


combinedIndexRemoved = [];
TFIndex = 1;
for myIndex=1:length(combinedIndex)
    if combinedIndex(myIndex)
        if TFRemoved(TFIndex)
            combinedIndexRemoved = [combinedIndexRemoved; true];
        else
            combinedIndexRemoved = [combinedIndexRemoved; false];
        end
        TFIndex = TFIndex + 1;
    else
        combinedIndexRemoved = [combinedIndexRemoved; false];
    end
end
combinedIndex = logical(combinedIndexRemoved);

for i=1:length(combinedIndex)
    if combinedIndex(i) == 0
        continue
    end
    if data.Psi(i) > 1000*median(data.Psi)
        combinedIndex(i) = 0;
    end
end
AllCount1HC = find(data.CellType(combinedIndex) == uniques(1,:));
AllCount1HGPS = find(data.CellType(combinedIndex) == uniques(2,:));
AllCount1Mild = find(data.CellType(combinedIndex) == uniques(3,:));
AllCount1Premanifest = find(data.CellType(combinedIndex) == uniques(4,:));
AllCount1Severe = find(data.CellType(combinedIndex) == uniques(5,:));
%
feature_vectors = [data.Dp(combinedIndex), data.Dtot(combinedIndex), data.Pp(combinedIndex), data.Psi(combinedIndex),...
    data.Pnp(combinedIndex), data.Dnp(combinedIndex), data.MSD10(combinedIndex), ...
    data.Sp(combinedIndex), data.Snp(combinedIndex)];
feature_names = [{'Dp'}, {'Dtot'},{'Pp'}, {'Psi'}, {'Pnp'}, {'Dnp'},...
    {'MSD10'}, {'Sp'}, {'Snp'}];
%% Prediction Intervals
temp_cellType = DataFileWell.CellType(combinedIndex);
patients = DataFileWell.CAPPatient(combinedIndex);
normalized_features = log(feature_vectors);
z = zscore(normalized_features);

featuresToUse = z;

featuresHC = featuresToUse(temp_cellType == "'HC'", :);
patientsHC = patients(temp_cellType == "'HC'");
patientsHCUnique = unique(patientsHC);
PfeaturesHC = zeros(size(patientsHCUnique, 1),9);
for i = 1:size(patientsHCUnique, 1)
   meanPatient = mean(featuresHC(patientsHC == patientsHCUnique(i), :));
   PfeaturesHC(i,:)= meanPatient;
end
featuresSevere = featuresToUse(temp_cellType == "'Severe'", :);
patientsSevere = patients(temp_cellType == "'Severe'");
patientsSevereUnique = unique(patientsSevere);
PfeaturesSevere = zeros(size(patientsSevereUnique, 1),9);
for i = 1:size(patientsSevereUnique, 1)
   meanPatient = mean(featuresSevere(patientsSevere == patientsSevereUnique(i), :));
   PfeaturesSevere(i,:)= meanPatient;
end
featuresMild = featuresToUse(temp_cellType == "'Mild'", :);
patientsMild = patients(temp_cellType == "'Mild'");
patientsMildUnique = unique(patientsMild);
PfeaturesMild = zeros(size(patientsMildUnique, 1),9);
for i = 1:size(patientsMildUnique, 1)
   meanPatient = mean(featuresMild(patientsMild == patientsMildUnique(i), :));
   PfeaturesMild(i,:)= meanPatient;
end
featuresPre = featuresToUse(temp_cellType == "'Premanifest'", :);
patientsPre = patients(temp_cellType == "'Premanifest'");
patientsPreUnique = unique(patientsPre);
PfeaturesPre = zeros(size(patientsPreUnique, 1),9);
for i = 1:size(patientsPreUnique, 1)
   meanPatient = mean(featuresPre(patientsPre == patientsPreUnique(i), :));
   PfeaturesPre(i,:)= meanPatient;
end
featuresHGPS = featuresToUse(temp_cellType == "'HGPS'", :);
patientsHGPS = patients(temp_cellType == "'HGPS'");
patientsHGPSUnique = unique(patientsHGPS);
PfeaturesHGPS = zeros(size(patientsHGPSUnique, 1),9);
for i = 1:size(patientsHGPSUnique, 1)
   meanPatient = mean(featuresHGPS(patientsHGPS == patientsHGPSUnique(i), :));
   PfeaturesHGPS(i,:)= meanPatient;
end

newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
confidenceZ = 1.5;

for feature=1:size(featuresToUse, 2)
    figure;
    hold on;
    
    % create a cell array to store features of each group
    featuresCell = {PfeaturesHC, PfeaturesPre, PfeaturesMild, PfeaturesSevere, PfeaturesHGPS};
    
    % create a cell array to store group names
    groupName = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};
    h = cell(1,5);

    for i = 1:5 % for each group
        disp("**********************");
        disp(groupName{i});
        
        featureUse = featuresCell{i};
        
        disp("Feature is:");
        disp(string(feature_names(feature)));
        
        currentFeatureMean = mean(featureUse(:, feature));
        currentFeatureSTD = std(featureUse(:, feature));
        
        disp("Lower Bound is:");
        lowerBound = currentFeatureMean - confidenceZ*currentFeatureSTD;
        disp(string(lowerBound));
        
        disp("Upper Bound is:");
        upperBound = currentFeatureMean + confidenceZ*currentFeatureSTD;
        disp(string(upperBound));
        
        % plot the histogram (you may want to adjust the number of bins)
        [N,edges] = histcounts(featureUse(:, feature), 'Normalization', 'probability');
        edges = edges(1:end-1) + diff(edges)/2;
%         h{i} = plot(edges, N, 'Color', newcolorsRGB{i});
        
        % plot the lower and upper bounds
        h{i} = line([lowerBound lowerBound], get(gca, 'YLim'), 'Color', newcolorsRGB{i}, 'LineStyle', '--');
        line([upperBound upperBound], get(gca, 'YLim'), 'Color', newcolorsRGB{i}, 'LineStyle', '--');
    end
    
    title(string(feature_names(feature)))
    xlabel('Feature Value')
    ylabel('Probability')
    legend([h{1} h{2} h{3} h{4} h{5}], groupName)
    
    hold off;
end

%% Z-Factor
temp_cellType = DataFileWell.CellType(combinedIndex);
normalized_features = log(feature_vectors);
z = zscore(normalized_features);

featuresToUse = z;


% featuresHC = featuresToUse(temp_cellType == "'HC'", :);
% featuresSevere = featuresToUse(temp_cellType == "'Severe'", :);
% featuresMild = featuresToUse(temp_cellType == "'Mild'", :);
% featuresPre = featuresToUse(temp_cellType == "'Premanifest'", :);
% featuresHGPS = featuresToUse(temp_cellType == "'HGPS'", :);

for feature=1:size(featuresToUse, 2)
    
    % create a cell array to store features of each group
    featuresCell = {PfeaturesHC, PfeaturesPre, PfeaturesMild, PfeaturesSevere, PfeaturesHGPS};
    
    % create a cell array to store group names
    groupName = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};
    h = cell(1,5);

    for i = 1:5 % for each group
        disp("**********************");
        disp("Group Outside is:");
        disp(groupName{i});
        
        featureUse = featuresCell{i};
        
        disp("Feature is:");
        disp(string(feature_names(feature)));
        
        currentFeatureMean = mean(featureUse(:, feature));
        currentFeatureSTD = std(featureUse(:, feature));
        for j = 1:5 % for each group
            if j == i
               continue; 
            end
            featureUseInside = featuresCell{j};

            disp("Group Inside is:");
            disp(groupName{j});

            insideFeatureMean = mean(featureUseInside(:, feature));
            insideFeatureSTD = std(featureUseInside(:, feature));
            
%             SeperationBand = (insideFeatureMean - 3*insideFeatureSTD) - (currentFeatureMean + 3*currentFeatureSTD);
            SeperationBand = 3*(insideFeatureSTD + currentFeatureSTD);
            DynamicRange = abs(insideFeatureMean - currentFeatureMean);
            zFactor = 1 - SeperationBand/DynamicRange;
            
            if zFactor > 0
                disp("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            end
            ("Z-Factor is:");
            disp(zFactor);
        end
    end
end
%%
tempType = data.CellType(combinedIndex);
tempPatient = data.CAPPatient(combinedIndex);
tempPatient(contains(string(tempPatient), "pro")) = "999'HGPS";
tempPatient(contains(string(tempPatient), "NA")) = "0'HC";

uniquePatients = unique(tempPatient);
uniquePatients = natsort(uniquePatients);

XLabels = str2double(regexp(string(uniquePatients),'\d*','Match', 'once'));
XLabels(XLabels==0) = 40;
XLabels(XLabels==999) = 161;
% if uniquePatients(end) == "'HC'"
%     uniquePatients = circshift(uniquePatients, 1);
% end
%
tileFirst = 1 + floor(size(feature_vectors,2) / 2);
meanListMatrix = zeros(size(feature_vectors,2), length(XLabels));
meanListMatrixGroup = zeros(size(feature_vectors,2), 5);

normalized_features = log(feature_vectors);
z = zscore(normalized_features);


for iBig=1:4:size(feature_vectors,2)
starts = iBig;
ends = iBig+3;
if ends > size(feature_vectors,2)
    ends = size(feature_vectors,2); 
end
fig = figure('Position', get(0, 'Screensize'));
t = tiledlayout(2,2, "TileSpacing", "loose");
for f=starts:ends

    meanList = [];
    meanStd = [];
    meanListZScore = [];
    meanListZScoreGroup = zeros(1,5);
    for i=1:length(uniquePatients)
        currentData = find(tempPatient == uniquePatients(i));
        currentFeatures = feature_vectors(currentData,f);
        currentFeaturesZScore = z(currentData,f);
        currentType = unique(tempType(currentData));
        if length(currentType) > 1
            disp('ERROR');
        end
        meanList(1+end) = mean(currentFeatures);
        meanListZScore(1+end) = mean(currentFeaturesZScore);
        if i == 1
            meanListZScoreGroup(1) = meanListZScoreGroup(1) + mean(currentFeaturesZScore);
        elseif i<8
                meanListZScoreGroup(2) = meanListZScoreGroup(2) + mean(currentFeaturesZScore);
        elseif i<17
                meanListZScoreGroup(3) = meanListZScoreGroup(3) + mean(currentFeaturesZScore);
        elseif i<22
                meanListZScoreGroup(4) = meanListZScoreGroup(4) + mean(currentFeaturesZScore);
        elseif i == 22
                meanListZScoreGroup(5) = meanListZScoreGroup(5) + mean(currentFeaturesZScore);
        end
        meanStd(1+end) = std(currentFeatures);
    end
    nexttile;
    currX = 1:length(meanList);
    XLabels = str2double(regexp(string(uniquePatients),'\d*','Match', 'once'));
    XLabels(XLabels==0) = 40;
    XLabels(XLabels==999) = 161;
    currX = XLabels;
    errorbar(currX, meanList,meanStd, 'bp', 'MarkerEdgeColor',[1 0 0],...
                'MarkerFaceColor',[0 .5 .5],...
                'LineWidth',0.25);
    meanListMatrix(f, :) = meanListZScore;
    meanListMatrixGroup(f,1) = meanListZScoreGroup(1);
    meanListMatrixGroup(f,2) = meanListZScoreGroup(2)/6;
    meanListMatrixGroup(f,3) = meanListZScoreGroup(3)/9;
    meanListMatrixGroup(f,4) = meanListZScoreGroup(4)/5;
    meanListMatrixGroup(f,5) = meanListZScoreGroup(5);

    hold on
    CAPS = unique(round(data.CAP(combinedIndex)));
    CAPS = CAPS(2:end-1)';
%     XT = [40, 60, 76, 94, 110, 119, 141, 161]; 
    XT = [40, CAPS, 161];
%     XLabels = ['HC', string(60), string(76), string(94), string(110), string(119), string(141), 'HGPS']; 
    XLabels = ['HC', string(CAPS), 'HGPS']; 
    set(gca, 'XTick', XT, 'XTickLabel', XLabels, 'fontweight','bold');
    pearson = corrcoef(currX,meanList);
    B = [currX(:) ones(size(currX(:)))] \ meanList(:);
    yfit = [currX(:) ones(size(currX(:)))]  * B;
    plot(currX, yfit, '-m')
    hold on
    fitresult = fit(currX,meanList','exp1');
    p11 = predint(fitresult,currX,0.8,'observation','off');
    hold on;
    plot(currX,p11,'k--')
    hold off
    xlim([min(XT)*0.8 max(XT)*1.05])
    ylim([-0.1 1.15*max(meanStd + meanList)])
    xlabel('CAP Score','fontweight','bold');
    ylabel(string(feature_names(f)),'fontweight','bold');
    if length(pearson) == 1
        title(strcat('Pearson Coefficient is', {' '}, num2str(pearson)),'fontweight','bold');
    else
        title(strcat('Pearson Coefficient is', {' '}, num2str(pearson(2))),'fontweight','bold');
    end

end
F    = getframe(fig);
exportgraphics(t,strcat(folderCurrent, '\_z_CAP_Correlations_', num2str(iBig), '.png'),'Resolution',500)
end

%% Velocity 
idxComb = combinedIndex;
xysM = xys_TRJs_with_HGPS(idxComb);
xys_full = xys_full_with_HGPS;

xysF = cell(length(xys_full),1);
for i=1:length(xys_full)
    current = xys_full{i};
    xysF{i} = xys_full{i}(:,[3,4]);
end
xysF = xysF(idxComb);

usedXY = xysM;
dt = 10 / 60; % time difference in hours
avgVelocity = zeros(size(usedXY, 1), 1);
avgDistance = zeros(size(usedXY, 1), 1);
meanderingIndex = zeros(size(usedXY, 1), 1);
outreachRatio = zeros(size(usedXY, 1), 1);
turningAngles = cell(size(xysF, 1), 1);

for i = 1:size(usedXY, 1)
    xy = usedXY{i, 1}; % extract the  matrix
    velocity = zeros(size(xy, 1)-1, 1); % preallocate for speed
    angles = zeros(size(xy, 1)-2, 1); % preallocate for speed
    distance = zeros(size(xy, 1)-1, 1); % preallocate for speed

    for j = 1:(size(xy, 1)-1)
        % calculate velocity between xy(j, :) and xy(j+1, :)
        dx = xy(j+1, 1) - xy(j, 1);
        dy = xy(j+1, 2) - xy(j, 2);
        distance(j) = sqrt(dx^2 + dy^2);
        velocity(j) = sqrt(dx^2 + dy^2) / dt;
    end
    avgDistance(i) = mean(distance);

    % store the average velocity for this row
    avgVelocity(i) = mean(velocity);
    
        % total path length
    pathLength = sum(sqrt(diff(xy(:, 1)).^2 + diff(xy(:, 2)).^2));

    % total displacement
    displacement = sqrt((xy(end, 1) - xy(1, 1))^2 + (xy(end, 2) - xy(1, 2))^2);

    % meandering index
    meanderingIndex(i) = displacement / pathLength;

    % maximum distance from the start point
    maxDist = max(sqrt((xy(:, 1) - xy(1, 1)).^2 + (xy(:, 2) - xy(1, 2)).^2));

    % outreach ratio
    outreachRatio(i) = maxDist / pathLength;
    
    for j = 1:(size(xy, 1)-2)
        % calculate the direction of movement at time j and j+1
        dir1 = atan2(xy(j+1, 2) - xy(j, 2), xy(j+1, 1) - xy(j, 1));
        dir2 = atan2(xy(j+2, 2) - xy(j+1, 2), xy(j+2, 1) - xy(j+1, 1));

        % calculate the turning angle
        angle = dir2 - dir1;

        % adjust the angle to the range -pi to pi
        angle = mod(angle + pi, 2*pi) - pi;

        % store the turning angle
        angles(j) = angle;
    end

    % store the turning angles for this row
    turningAngles{i} = angles;
end

weights = [1/3, 1/3, 1/3]; % equally important
% normalize to 0-1
normalizedVelocity = (avgVelocity - min(avgVelocity)) / (max(avgVelocity) - min(avgVelocity));
normalizedMeanderingIndex = (meanderingIndex - min(meanderingIndex)) / (max(meanderingIndex) - min(meanderingIndex));
normalizedOutreachRatio = (outreachRatio - min(outreachRatio)) / (max(outreachRatio) - min(outreachRatio));
% combine into a single score
combinedScore = weights(1) * normalizedVelocity + weights(2) * normalizedMeanderingIndex + weights(3) * normalizedOutreachRatio;

meanTurningAngle = cellfun(@mean, turningAngles);
medianTurningAngle = cellfun(@median, turningAngles);
stdTurningAngle = cellfun(@std, turningAngles);
circVarTurningAngle = cellfun(@(x) 1 - abs(sum(exp(1i * x))) / length(x), turningAngles);
weightsCirc = [0.15, 0.85]; % equally important
normalizedCircVarTurningAngle = (circVarTurningAngle - min(circVarTurningAngle)) / (max(circVarTurningAngle) - min(circVarTurningAngle));
combinedScoreCirc = weights(1) * normalizedVelocity + weights(2) * normalizedCircVarTurningAngle;

whatToPlot = circVarTurningAngle;

meanList = [];
meanStd = [];
meanListZScore = [];
meanListZScoreGroup = zeros(1,5);
for i=1:length(uniquePatients)
    currentData = find(tempPatient == uniquePatients(i));
    currentFeatures = whatToPlot(currentData);
    currentType = unique(tempType(currentData));
    if length(currentType) > 1
        disp('ERROR');
    end
    meanList(1+end) = mean(currentFeatures);
    meanStd(1+end) = std(currentFeatures);
end
nexttile;
currX = 1:length(meanList);
XLabels = str2double(regexp(string(uniquePatients),'\d*','Match', 'once'));
XLabels(XLabels==0) = 40;
XLabels(XLabels==999) = 161;
currX = XLabels;
errorbar(currX, meanList,meanStd, 'bp', 'MarkerEdgeColor',[1 0 0],...
            'MarkerFaceColor',[0 .5 .5],...
            'LineWidth',0.25);

hold on
CAPS = unique(round(data.CAP(combinedIndex)));
CAPS = CAPS(2:end-1)';
%     XT = [40, 60, 76, 94, 110, 119, 141, 161]; 
XT = [40, CAPS, 161];
%     XLabels = ['HC', string(60), string(76), string(94), string(110), string(119), string(141), 'HGPS']; 
XLabels = ['HC', string(CAPS), 'HGPS']; 
set(gca, 'XTick', XT, 'XTickLabel', XLabels, 'fontweight','bold');
pearson = corrcoef(currX,meanList);

% B = [currX(:) ones(size(currX(:)))] \ meanList(:);
% yfit = [currX(:) ones(size(currX(:)))]  * B;

p = polyfit(currX', meanList', 2);
yfit = polyval(p, currX);

plot(currX, yfit, '-m')
hold on
fitresult = fit(currX,meanList','exp1');
p11 = predint(fitresult,currX,0.8,'observation','off');
hold on;
plot(currX,p11,'k--')
hold off
xlim([min(XT)*0.8 max(XT)*1.05])
ylim([-0.1 1.15*max(meanStd + meanList)])
xlabel('CAP Score','fontweight','bold');
ylabel('Score','fontweight','bold');
if length(pearson) == 1
    title(strcat('Pearson Coefficient is', {' '}, num2str(pearson)),'fontweight','bold');
else
    title(strcat('Pearson Coefficient is', {' '}, num2str(pearson(2))),'fontweight','bold');
end

%% Plot the confusion matrix
figure; % Create a new figure
imagesc(1:5, 1:size(feature_vectors,2), meanListMatrixGroup); % Plot the matrix
colorbar; % Add a colorbar to show the intensity scale
xlabel('Group');
ylabel('Feature');
title('Correlations');
ax = gca;
ax.XTick = 1:5;
ax.YTick = 1:size(feature_vectors,2);
xLabels = {'HC',"Premanifest","Mild","Severe","HGPS"};
ax.XTickLabel = string(xLabels);
ax.YTickLabel = string(feature_names);
% Define the color gradient
numColors = 256; % The number of colors in the colormap
% Define the most red (#DB4437), white (#FFFFFF), and most green (#0F9D58) colors in normalized RGB format
mostRed = [219, 68, 55] / 255;
white = [255, 255, 255] / 255;
mostGreen = [15, 157, 88] / 255;
lessGreen = [210, 252, 232] / 255;
leesRed = [252, 238, 237] / 255;
% Define the control points for the color gradient (positions between 0 and 1)
% Adjust the control points to create a sharper transition around zero
controlPoints = [0, 0.45, 0.5, 0.55, 1];
colors = [mostRed; leesRed; white; lessGreen; mostGreen];
% Create a non-linear gradient using interp1
xq = linspace(0, 1, numColors);
customColormap = interp1(controlPoints, colors, xq, 'pchip'); % 'pchip' method provides smooth transitions
colormap(customColormap);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\GroupCorrelationsHeatmap.png'))


% Plot the confusion matrix
figure; % Create a new figure
imagesc(1:22, 1:size(feature_vectors,2), meanListMatrix); % Plot the matrix
colorbar; % Add a colorbar to show the intensity scale
xlabel('CAP Score');
ylabel('Feature');
title('Correlations');
ax = gca;
ax.XTick = 1:22;
ax.YTick = 1:size(feature_vectors,2);
xLabels = {'HC',"60",    "69",    "75",    "76a","76b",    "84",    "91",...
        "94",    "96",    "100",    "104",    "106",    "110",    "111a", "111b",   "115",    "119",...
    "126",    "129",    "141",    "HGPS"};
ax.XTickLabel = string(xLabels);
ax.YTickLabel = string(feature_names);
% Define the color gradient
numColors = 256; % The number of colors in the colormap
% Define the most red (#DB4437), white (#FFFFFF), and most green (#0F9D58) colors in normalized RGB format
mostRed = [219, 68, 55] / 255;
white = [255, 255, 255] / 255;
mostGreen = [15, 157, 88] / 255;
lessGreen = [210, 252, 232] / 255;
leesRed = [252, 238, 237] / 255;
% Define the control points for the color gradient (positions between 0 and 1)
% Adjust the control points to create a sharper transition around zero
controlPoints = [0, 0.45, 0.5, 0.55, 1];
colors = [mostRed; leesRed; white; lessGreen; mostGreen];
% Create a non-linear gradient using interp1
xq = linspace(0, 1, numColors);
customColormap = interp1(controlPoints, colors, xq, 'pchip'); % 'pchip' method provides smooth transitions
colormap(customColormap);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\CAPCorrelationsHeatmap.png'))

useMe = z;
% Compute the correlation matrix
correlationMatrix = corrcoef(useMe);
% Create a heatmap of the correlation matrix
figure;
imagesc(correlationMatrix);
colorbar; % Add a colorbar to show the intensity scale
colormap(jet);
% Set the axis labels
xlabel('Feature');
ylabel('Feature');
title('Correlation Matrix Heatmap');
% Set the x and y axis values to match your feature names
ax = gca;
ax.XTickLabel = string(feature_names);
ax.YTickLabel = string(feature_names);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\FeaturCorrelations.png'))

% Compute the cross-correlation matrix
numFeatures = size(useMe, 2);
crossCorrMatrix = zeros(numFeatures);
for i = 1:numFeatures
    for j = 1:numFeatures
        tempCrossCorr = xcorr(useMe(:, i), useMe(:, j), 'coeff'); % 'coeff' normalizes the cross-correlation
        crossCorrMatrix(i, j) = max(abs(tempCrossCorr)); % Store the maximum magnitude of cross-correlation
    end
end
% Create a heatmap of the correlation matrix
figure;
imagesc(crossCorrMatrix);
colorbar; % Add a colorbar to show the intensity scale
colormap(jet);
% Set the axis labels
xlabel('Feature');
ylabel('Feature');
title('Cross-Correlation Matrix Heatmap');
% Set the x and y axis values to match your feature names
ax = gca;
ax.XTickLabel = string(feature_names);
ax.YTickLabel = string(feature_names);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\FeaturCrossCorrelations.png'))

% close all
%
% xysMain=get_trajfile; %(DATA.mat)
% xysMain_progeria=get_trajfile; %(DATA.xlsx)

%%%
idxComb = combinedIndex;
xysMain = xys_TRJs_with_HGPS;
xys_full = xys_full_with_HGPS;
% xys = xysMain(idxComb);

newcolors = {'#4285F4','#0F9D58','#F4B400','#DB4437' '#330072', '#016773'};
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};


xys = cell(length(xys_full),1);
for i=1:length(xys_full)
    current = xys_full{i};
    xys{i} = xys_full{i}(:,[3,4]);
end
xys = xys(idxComb);

% MSD vs TLAG
% close all force
param.showfig=1;
param.saveres=0;
param.markertype='.-';
param.outfigurenum=301;
param.dim=2;
param.linear=0;
param.aviv=0;

temp_cellType = DataFileWell.CellType(idxComb);

a1Pol = xys(temp_cellType == "'HC'");
a2Pol = xys(temp_cellType == "'Severe'");
a3Pol = xys(temp_cellType == "'Mild'");
a4Pol = xys(temp_cellType == "'Premanifest'");
a5Pol = xys(temp_cellType == "'HGPS'");

get_MSD(a1Pol, 10, param);
hold on;
param.markertype='.-';
get_MSD(a4Pol, 10, param);
hold on;
param.markertype='.-';
get_MSD(a3Pol, 10, param);
hold on;
param.markertype='.-'; 
get_MSD(a2Pol, 10, param);
hold on;
param.markertype='.-'; 
get_MSD(a5Pol, 10, param);
hold on;

colororder(newcolors)
param.markertype='k:';
param.linear=1;
get_MSD(a1Pol, 10, param);
hold off;
title('Mean Squared Displacement');

ylabel('MSD (µm2)') ;
xlabel('Time lag (min)');

legend('HC','Premanifest' ,'Mild', 'Severe','HGPS','Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_MSD.png'))
%
%%% MSD vs TLAG PATIENT
close all force

param.showfig=1;
param.saveres=0;
param.aviv=1;
param.outfigurenum=301;
param.dim=2;
param.linear=0;

temp_Patient = DataFileWell.CAP(idxComb);
cidk = jet(length(unique(temp_Patient)));
name = sort(unique(temp_Patient), 'descend');
for i=1:length(name)
    
    temp_pat = i;

    clust_idx = find(temp_Patient == name(i));
    check = find(temp_Patient == name(temp_pat));
 
    aPol = xys(check);
    
    param.aviv2 = 0;
    param.markertype=cidk(length(unique(temp_Patient)) + 1 - i,:);
    get_MSD(aPol, 10, param);
    hold on;

end



hold on;
param.markertype='k:';
param.linear=1;
get_MSD(aPol, 10, param);

title('Mean Squared Displacement by CAP Scores');
ylabel('MSD (µm2)') ;
xlabel('Time lag (min)');

name = string(name);
name(1) = 'HGPS';
name(end) = 'HC';
legend(name, 'Location', 'bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_MSDPatients.png'))
%
% dR Polarity
close all
param.showfig=1;
param.saveres=0;
param.dim=2;
param.outfigurenum=305;
param.binnum=35;

param.markertype=newcolorsRGB{1};
get_dR_polarity(a1Pol, 1, param);
hold on;
param.markertype=newcolorsRGB{2};
get_dR_polarity(a4Pol, 1, param);
hold on;
param.markertype=newcolorsRGB{3};
get_dR_polarity(a3Pol, 1, param);
hold on;
param.markertype=newcolorsRGB{4};
get_dR_polarity(a2Pol, 1, param);
hold on;
param.markertype=newcolorsRGB{5};
get_dR_polarity(a5Pol, 1, param);
hold on;
% colororder(newcolors)

title('dR Polarity');

g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_Theta.png'))
% legend('HC','Premanifest' ,'Mild', 'Severe','Progeria','Location','bestoutside');
%
% dR Polarity Patient
close all
param.showfig=1;
param.saveres=0;
param.dim=2;
param.outfigurenum=305;
param.binnum=35;


temp_Patient = DataFileWell.CAP(idxComb);
cidk = jet(length(unique(temp_Patient)));
name = sort(unique(temp_Patient), 'descend');
for i=1:length(name)
    
    temp_pat = i;
    clust_idx = find(temp_Patient == name(i));
    check = find(temp_Patient == name(temp_pat));
 
    aPol = xys(check);

    param.markertype=cidk(length(unique(temp_Patient)) + 1 - i,:);
    get_dR_polarity(aPol, 1, param);
    hold on;
    
    
end

title('dR Polarity');
name = string(name);
name(1) = 'HGPS';
name(end) = 'HC';
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_ThetaPatients.png'))
% legend(name, 'Location', 'bestoutside');
% legend('HC','Premanifest' ,'Mild', 'Severe','Location','bestoutside');
%
%% dTheta PDF
close all
param.showfig=1;
param.saveres=0;
param.dim=2;
param.outfigurenum=405;
param.binnum=6;
param.markertype='-';
tloi = [1 2 3 4 5 6];
cid=copper(length(tloi));
% param.markertype=newcolorsRGB{1};
% figure();
get_dtheta_PDF(a1Pol, tloi, param);
% legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
title('HC dθ PDF');
xlabel('dθ') ;
ylabel('Occurrence');
    axis square;
    box on
        colorMatrix = repmat(linspace(1, length(tloi), 100), length(tloi), 1);
        color_ax = axes('Position', [0.7, 0.1, 0.2, 0.05]);
        contourf(color_ax, colorMatrix, length(tloi), 'LineStyle', 'none');
        colormap(color_ax, cid);
        axis off;
        color_bar_title = title(color_ax, '0 - 60 Minutes', 'Units', 'normalized', 'Position', [0.5, 1.3, 0],'FontSize', 18, 'FontWeight', 'bold');

g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_Theta_PDF_HC.png'))
% hold on;
% param.markertype=newcolorsRGB{2};
% figure();
param.outfigurenum = param.outfigurenum + 1;
get_dtheta_PDF(a4Pol, tloi, param);
% legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
title('Premanifest dθ PDF');
xlabel('dθ') ;
ylabel('Occurrence');
    axis square;
    box on
        colorMatrix = repmat(linspace(1, length(tloi), 100), length(tloi), 1);
        color_ax = axes('Position', [0.7, 0.1, 0.2, 0.05]);
        contourf(color_ax, colorMatrix, length(tloi), 'LineStyle', 'none');
        colormap(color_ax, cid);
        axis off;
        color_bar_title = title(color_ax, '0 - 60 Minutes', 'Units', 'normalized', 'Position', [0.5, 1.3, 0],'FontSize', 18, 'FontWeight', 'bold');

        g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_Theta_PDF_Premanifest.png'))
% hold on;
% param.markertype=newcolorsRGB{3};
% figure();
param.outfigurenum = param.outfigurenum + 1;
get_dtheta_PDF(a3Pol, tloi, param);
% legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
title('Mild dθ PDF');
xlabel('dθ') ;
ylabel('Occurrence');
    axis square;
    box on
        colorMatrix = repmat(linspace(1, length(tloi), 100), length(tloi), 1);
        color_ax = axes('Position', [0.7, 0.1, 0.2, 0.05]);
        contourf(color_ax, colorMatrix, length(tloi), 'LineStyle', 'none');
        colormap(color_ax, cid);
        axis off;
        color_bar_title = title(color_ax, '0 - 60 Minutes', 'Units', 'normalized', 'Position', [0.5, 1.3, 0],'FontSize', 18, 'FontWeight', 'bold');

        g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_Theta_PDF_Mild.png'))
% hold on;
% param.markertype=newcolorsRGB{4};
% figure();
param.outfigurenum = param.outfigurenum + 1;
get_dtheta_PDF(a2Pol, tloi, param);
% legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
title('Severe dθ PDF');
xlabel('dθ') ;
ylabel('Occurrence');
    axis square;
    box on
        colorMatrix = repmat(linspace(1, length(tloi), 100), length(tloi), 1);
        color_ax = axes('Position', [0.7, 0.1, 0.2, 0.05]);
        contourf(color_ax, colorMatrix, length(tloi), 'LineStyle', 'none');
        colormap(color_ax, cid);
        axis off;
        color_bar_title = title(color_ax, '0 - 60 Minutes', 'Units', 'normalized', 'Position', [0.5, 1.3, 0],'FontSize', 18, 'FontWeight', 'bold');
        g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_Theta_PDF_Severe.png'))
% hold on;
% param.markertype=newcolorsRGB{5};
% figure();
param.outfigurenum = param.outfigurenum + 1;
get_dtheta_PDF(a5Pol, tloi, param);
% legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
title('HGPS dθ PDF');
xlabel('dθ') ;
ylabel('Occurrence');
    axis square;
    box on
        colorMatrix = repmat(linspace(1, length(tloi), 100), length(tloi), 1);
        color_ax = axes('Position', [0.7, 0.1, 0.2, 0.05]);
        contourf(color_ax, colorMatrix, length(tloi), 'LineStyle', 'none');
        colormap(color_ax, cid);
        axis off;
        color_bar_title = title(color_ax, '0 - 60 Minutes', 'Units', 'normalized', 'Position', [0.5, 1.3, 0],'FontSize', 18, 'FontWeight', 'bold');

        g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_Theta_PDF_HGPS.png'))
% hold on;
% colororder(newcolors)


% legend('HC','Premanifest' ,'Mild', 'Severe','Progeria','Location','bestoutside');
%
% dTheta PDF Patient
% close all
% param.showfig=1;
% param.saveres=0;
% param.dim=2;
% param.outfigurenum=305;
% param.binnum=35;
% 
% 
% temp_Patient = DataFileWell.CAP(idxComb);
% cidk = jet(length(unique(temp_Patient)));
% name = sort(unique(temp_Patient), 'descend');
% for i=1:length(name)
%     
%     temp_pat = i;
%     clust_idx = find(temp_Patient == name(i));
%     check = find(temp_Patient == name(temp_pat));
%  
%     aPol = xys(check);
% 
%     param.markertype=cidk(length(unique(temp_Patient)) + 1 - i,:);
%     get_dtheta_PDF(aPol, tloi, param);
%     hold on;
%     
%     
% end
% 
% title('dθ PDF');
% name = string(name);
% name(1) = 'HGPS';
% name(end) = 'HC';
%     axis square;
%     box on
% g = gcf;
% g.WindowState = 'maximized';
% saveas(g, strcat(folderCurrent, '\z_ThetaPatients_PDF.png'))
% legend(name, 'Location', 'bestoutside');
% legend('HC','Premanifest' ,'Mild', 'Severe','Location','bestoutside');
%% dTheta PDF by Patient

close all
param.showfig=1;
param.saveres=0;
param.dim=2;
param.outfigurenum=805;
param.binnum=6;
param.markertype='-';
tloi = [1 2 3 4 5 6];
cid=copper(length(tloi));

temp_Patient = DataFileWell.CAP(idxComb);
name = sort(unique(temp_Patient), 'descend');
    
    namestr = string(name);
    namestr(1) = 'HGPS';
    namestr(end) = 'HC';
    
for i=1:length(unique(temp_Patient))
    
    temp_pat = i;

    clust_idx = find(temp_Patient == name(i));
    check = find(temp_Patient == name(temp_pat));
 
    aPol = xys(check);
    
    param.outfigurenum = param.outfigurenum + 1;
    get_dtheta_PDF(aPol, tloi, param);
    
    if namestr(i) == 'HC'
        title(strcat('dθ PDF of ', {' '}, namestr(i)));
    elseif namestr(i) == 'HGPS'
        title(strcat('dθ PDF of ', {' '}, namestr(i)));
    else
        title(strcat('dθ PDF of CAP ', {' '}, namestr(i)));
    end
    xlabel('dθ') ;
    ylabel('Occurrence');
    axis square;
    box on
    colorMatrix = repmat(linspace(1, length(tloi), 100), length(tloi), 1);
    color_ax = axes('Position', [0.7, 0.1, 0.2, 0.05]);
    contourf(color_ax, colorMatrix, length(tloi), 'LineStyle', 'none');
    colormap(color_ax, cid);
    axis off;
    color_bar_title = title(color_ax, '0 - 60 Minutes', 'Units', 'normalized', 'Position', [0.5, 1.3, 0],'FontSize', 18, 'FontWeight', 'bold');

    g = gcf;
    g.WindowState = 'maximized';
    saveas(g, strcat(folderCurrent, '\z_Theta_PDF_Patients_', namestr(i),'.png'))
end

%%
% ACF vs TLAG
close all force
param.dim=2; 
param.tlag=1;
param.saveres=0;
param.showfig=1;       
param.markertype=''; % HC
param.outfigurenum=303;   
param.aviv=1;
param.aviv3 = 0;

xys = cell(length(xys_full),1);
for i=1:length(xys_full)
    current = xys_full{i};
    xys{i} = xys_full{i}(:,[3,4]);
end
xys = xys(idxComb);

a1Pol = xys(temp_cellType == "'HC'");
a2Pol = xys(temp_cellType == "'Severe'");
a3Pol = xys(temp_cellType == "'Mild'");
a4Pol = xys(temp_cellType == "'Premanifest'");
a5Pol = xys(temp_cellType == "'HGPS'");

get_ACF1(a1Pol, 10, param);
hold on;
get_ACF1(a4Pol, 10, param);
hold on;
get_ACF1(a3Pol, 10, param);
hold on;
get_ACF1(a2Pol, 10, param);
hold on;
get_ACF1(a5Pol, 10, param);
hold on;

% newcolors = {'#4285F4','#4285F4', '#0F9D58','#0F9D58','#F4B400','#F4B400','#DB4437','#DB4437', '#330072', '#330072','#016773', '#016773'};
% newcolorsRGB = {[66 133 244]/255,[0 0 0],[15 157 88]/255,[0 0 0],[244 180 0]/255,[0 0 0],[219 68 55]/255, [0 0 0],[51 0 114]/255, [0 0 0]};
newcolors = {'#4285F4','#0F9D58','#F4B400','#DB4437' '#330072', '#016773'};
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colororder(newcolors)

title('Auto-Correlation Function');
ylabel('ACF') ;
xlabel('Time lag (min)');
% legend('HC','line','Premanifest' ,'line','Mild', 'line','Severe','line','HGPS','line','Location','bestoutside');
legend('HC','Premanifest','Mild','Severe','HGPS','Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_ACF.png'))
%%
% ACF vs TLAG PATIENT
close all force
param.dim=2; 
param.tlag=1;
param.saveres=0;
param.showfig=1;       
param.markertype2='none'; % HC
param.outfigurenum=303;   
param.aviv=1;


temp_Patient = DataFileWell.CAP(idxComb);
cidk = jet(length(unique(temp_Patient)));
name = sort(unique(temp_Patient), 'descend');
for i=1:length(unique(temp_Patient))
    
    temp_pat = i;

    clust_idx = find(temp_Patient == name(i));
    check = find(temp_Patient == name(temp_pat));
 
    aPol = xys(check);
    
    param.aviv2 = 0;
    param.aviv3 = 1;
    param.markertype=cidk(length(unique(temp_Patient)) + 1 - i,:);
    get_ACF_Patients(aPol, 10, param);
    hold on;

    
end

title('Auto-Correlation Function by CAP Scores');
ylabel('ACF') ;
xlabel('Time lag (min)');

name = string(name);
name(1) = 'HGPS';
name(end) = 'HC';
legend(name, 'Location', 'bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_ACFPatients.png'))
%
%% dR PDF
close all force
param.dxmax=50;
param.binn=70;            
param.showfig=1;
param.saveres=0;
param.dim=2;            
param.outfigurenum=302;
param.markertype='o'; % HC
param.aviv=0;
param.aviv3 = 1;

param.markertype=':';
get_dR_PDF(a1Pol, 10, param);
hold on;
param.markertype=':'; % HD
get_dR_PDF(a4Pol, 10, param);
hold on;
param.markertype=':'; % HD
get_dR_PDF(a3Pol, 10, param);
hold on;
param.markertype=':'; % HD
get_dR_PDF(a2Pol, 10, param);
hold on;
param.markertype=':'; % HD
get_dR_PDF(a5Pol, 10, param);
hold on;
colororder(newcolors)

title('PDF Cellular Displacements');
ylabel('Occurrence') ;
xlabel('Displacement (µm)');
legend('HC','Premanifest' ,'Mild', 'Severe','HGPS','Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_PDF.png'))
%
%%% dR PDF PATIENT
close all force
param.dxmax=50;
param.binn=70;            
param.showfig=1;
param.saveres=0;
param.dim=2;            
param.outfigurenum=302;
param.markertype='bo'; % HC  
param.aviv=1;

temp_Patient = DataFileWell.CAP(idxComb);
cidk = jet(length(unique(temp_Patient)));
name = sort(unique(temp_Patient), 'descend');
for i=1:length(unique(temp_Patient))
    
    temp_pat = i;

    clust_idx = find(temp_Patient == name(i));
    check = find(temp_Patient == name(temp_pat));
 
    aPol = xys(check);
   
    param.aviv2=1;
    param.markertype=cidk(length(unique(temp_Patient)) + 1 - i,:);
    get_dR_PDF(aPol, 10, param);
    hold on;

end

title('PDF Cellular Displacements by CAP Scores');
ylabel('Occurrence') ;
xlabel('Displacement (µm)');

name = string(name);
name(1) = 'HGPS';
name(end) = 'HC';
legend(name, 'Location', 'bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_PDFPatients.png'))
%
%%% Draw TRJ

normalized_features = log(feature_vectors);
z = zscore(normalized_features);

xys = xysMain(idxComb);
a1Pol = xys(temp_cellType == "'HC'");
a2Pol = xys(temp_cellType == "'Severe'");
a3Pol = xys(temp_cellType == "'Mild'");
a4Pol = xys(temp_cellType == "'Premanifest'");
a5Pol = xys(temp_cellType == "'HGPS'");
%%% Show Random Cells
close all force
numOfRandoms = 50;

for SMR=1:5
    fig = figure('Position', get(0, 'Screensize'));
    t = tiledlayout(3,1, 'TileSpacing','tight');
    if SMR == 1
        cells = a2Pol;
        titleSMR = {"Trajectories", "Severe"};
        color = newcolors{4};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
        out1 = randperm(size(cells,1),numOfRandomsIn);
        out1All = randperm(size(cells,1),size(cells,1));
    elseif SMR == 2
        cells = a3Pol;
        titleSMR = {"Trajectories", "Mild"};
        color = newcolors{3};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
        out2 = randperm(size(cells,1),numOfRandomsIn);
        out2All = randperm(size(cells,1),size(cells,1));
    elseif SMR == 3
        cells = a4Pol;
        titleSMR = {"Trajectories", "Premanifest"};
        color = newcolors{2};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
        out3 = randperm(size(cells,1),numOfRandomsIn);
        out3All = randperm(size(cells,1),size(cells,1));
    elseif SMR == 4
        cells = a1Pol;
        titleSMR = {"Trajectories", "HC"};
        color = newcolors{1};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
        out4 = randperm(size(cells,1),numOfRandomsIn);
        out4All = randperm(size(cells,1),size(cells,1));
    elseif SMR == 5
        cells = a5Pol;
        titleSMR = {"Trajectories", "HGPS"};
        color = newcolors{5};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
        out5 = randperm(size(cells,1),numOfRandomsIn);
        out5All = randperm(size(cells,1),size(cells,1));
    end
    counter_x = 0;
    counter_y = 0;
    space = 500;
    nexttile
    hold on
    for j=1:numOfRandomsIn
        if SMR == 1
            currNum = out1(j);
        elseif SMR == 2
            currNum = out2(j);
        elseif SMR == 3
            currNum = out3(j);
        elseif SMR == 4
            currNum = out4(j);
        elseif SMR == 5
            currNum = out5(j);
        end
        % getting sepcific cell trajectory
        xy=cells{currNum}; 

        % zero cell position at first frame.
        xy=xy-ones(size(xy(:,1)))*xy(1,:);
        plot(xy(:,1)+counter_x*space,xy(:,2)+counter_y*space,'Color',color, 'LineWidth', 2);
        counter_x = counter_x + 1;
        if mod(counter_x,10)==0
            counter_x = 0;
            counter_y = counter_y + 1;
        end
    end
    set(gca,'xtick',[],'ytick',[]);
    title(titleSMR);
    axis square;
    box on
    hold off
    ylim([-space space*5])
    xlim([-space space*10])


% Show Sunplot of Cells


    if SMR == 1
        cells = a2Pol;
        titleSMR = strcat("Severe Sunplot Trajectories");
        color = newcolors{4};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 2
        cells = a3Pol;
        titleSMR = strcat("Mild Sunplot Trajectories");
        color = newcolors{3};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 3
        cells = a4Pol;
        titleSMR = strcat("Premanifest Sunplot Trajectories");
        color = newcolors{2};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 4
        cells = a1Pol;
        titleSMR = strcat("HC Sunplot Trajectories");
        color = newcolors{1};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 5
        cells = a5Pol;
        titleSMR = strcat("HGPS Sunplot Trajectories");
        color = newcolors{5};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    end
    nexttile
    cidk = jet(numOfRandomsIn);
    for k=1:1:numOfRandomsIn
        if SMR == 1
            currNum = out1(k);
        elseif SMR == 2
            currNum = out2(k);
        elseif SMR == 3
            currNum = out3(k);
        elseif SMR == 4
            currNum = out4(k);
        elseif SMR == 5
            currNum = out5(k);
        end
         xy=cells{currNum};         
         xy=xy-ones(size(xy(:,1)))*xy(1,:); % zero cell position at first frame.
         xy1 = smooth(xy(:,1));
         xy2 = smooth(xy(:,2));
         plot(xy1,xy2,'-','color',[cidk(k,:) 1],'linewidth',2); hold on; 
    end
%     title(titleSMR);
    axis square
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-350 350])
    xlim([-350 350])
    text(90, 300, strcat("N = ", string(length(cells))),'FontSize', 16, 'FontWeight', 'bold');

    
    % New Plt medged

    if SMR == 1
        cells = a2Pol;
        titleSMR = strcat("Severe Merged Trajectory");
        color = newcolors{4};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 2
        cells = a3Pol;
        titleSMR = strcat("Mild Merged Trajectory");
        color = newcolors{3};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 3
        cells = a4Pol;
        titleSMR = strcat("Premanifest Merged Trajectory");
        color = newcolors{2};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 4
        cells = a1Pol;
        titleSMR = strcat("HC Merged Trajectory");
        color = newcolors{1};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 5
        cells = a5Pol;
        titleSMR = strcat("HGPS Merged Trajectory");
        color = newcolors{5};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    end
    nexttile
    cidk = jet(numOfRandomsIn);
    xy = 0;
    for k=1:1:size(cells,1)
        if SMR == 1
            currNum = out1All(k);
        elseif SMR == 2
            currNum = out2All(k);
        elseif SMR == 3
            currNum = out3All(k);
        elseif SMR == 4
            currNum = out4All(k);
        elseif SMR == 5
            currNum = out5All(k);
        end
         xyC=cells{currNum};         
         xy=xy + xyC-ones(size(xyC(:,1)))*xyC(1,:); % zero cell position at first frame.
    end
    xy = xy / size(cells,1);
    xy1 = smooth(xy(:,1));
    xy2 = smooth(xy(:,2));
    plot(xy1,xy2,'-','color',color, 'linewidth',2);
%     title(titleSMR);
    axis square
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-15 15])
    xlim([-15 15])
    exportgraphics(t,strcat(folderCurrent, 'a_', string(SMR),'_TRJFirst.jpg'),'Resolution',500)

end
%
close all force
%%% Show Random Cells Per CAP
numOfRandoms = 50;

temp_Patient = DataFileWell.CAP(idxComb);
name = sort(unique(temp_Patient), 'descend');
for i=1:length(unique(temp_Patient))
    fig = figure('Position', get(0, 'Screensize'));
    t = tiledlayout(3,1, 'TileSpacing','tight');
    color = newcolors{4};
    if name(i) == 0
        color = newcolors{1};
    elseif name(i) < 90
        color = newcolors{2};
    elseif name(i) < 114
        color = newcolors{3};
    elseif name(i) == 999
        color = newcolors{5};
    end
    temp_pat = i;
    check = find(temp_Patient == name(temp_pat));
    aPol = xys(check);
    if size(aPol,1) < numOfRandoms
        numOfRandomsIn = size(aPol,1);
    else
        numOfRandomsIn = numOfRandoms;
    end
    out = randperm(size(aPol,1),numOfRandomsIn);
    outAll = randperm(size(aPol,1),size(aPol,1));
    counter_x = 0;
    counter_y = 0;
    space = 500;
    nexttile
    hold on
    for j=1:numOfRandomsIn
        currNum = out(j);
        % getting sepcific cell trajectory
        xy=aPol{currNum}; 

        % zero cell position at first frame.
        xy=xy-ones(size(xy(:,1)))*xy(1,:);
        plot(xy(:,1)+counter_x*space,xy(:,2)+counter_y*space,'Color',color, 'LineWidth', 2);
        counter_x = counter_x + 1;
        if mod(counter_x,10)==0
            counter_x = 0;
            counter_y = counter_y + 1;
        end
    end
    set(gca,'xtick',[],'ytick',[]);
    stringName = name(i);
    if stringName == 0
        title({"Trajectories", "HC"});
    elseif stringName == 999
        title({"Trajectories", "HGPS"});
    else
        title({"Trajectories of CAP:", string(round(name(i)*100)/100)});
    end
    axis square;
    box on
    hold off
    ylim([-space space*5])
    xlim([-space space*10])
% Show Sunplot of Cells
    nexttile
    cidk = jet(numOfRandomsIn);
    for k=1:1:numOfRandomsIn
         currNum = out(k);
         xy=aPol{currNum};         
         xy=xy-ones(size(xy(:,1)))*xy(1,:); % zero cell position at first frame.
           xy1 = smooth(xy(:,1));
         xy2 = smooth(xy(:,2));
         plot(xy1,xy2,'-','color',[cidk(k,:) 1],'linewidth',2); hold on; 
    end
%     if stringName == 0
%         title(strcat("Sunplot Trajectories for HC Cells"));
%     elseif stringName == 999
%         title(strcat("Sunplot Trajectories for HGPS Cells"));
%     else
%         title(strcat("Sunplot Trajectories for CAP Score: ", string(name(i))));
%     end
    axis square
        box on
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-350 350])
    xlim([-350 350])
    text(90, 300, strcat("N = ", string(length(aPol))),'FontSize', 16, 'FontWeight', 'bold');

    % Show Mrged of Cells
    nexttile
    cidk = jet(numOfRandomsIn);
    xy = 0;
    for k=1:1:numOfRandomsIn
         currNum = out(k);
         xyC=aPol{currNum};         
         xy=xy + xyC-ones(size(xyC(:,1)))*xyC(1,:); % zero cell position at first frame.
        
    end
    xy = xy / numOfRandomsIn;
    xy1 = smooth(xy(:,1));
    xy2 = smooth(xy(:,2));
    plot(xy1,xy2,'-','color',color, 'linewidth',2);
%     if stringName == 0
%         title(strcat("Merged Trajectory for HC Cells"));
%     elseif stringName == 999
%         title(strcat("Merged Trajectory for HGPS Cells"));
%     else
%         title(strcat("Merged Trajectory for CAP Score: ", string(name(i))));
%     end
    axis square
        box on
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-30 30])
    xlim([-30 30])
    exportgraphics(t,strcat(folderCurrent, 'b_', string(i),'_TRJCAP.jpg'),'Resolution',500)

end
%
% close all force
% PCA Features
features_T = z';
[coeff, score, latent, tsquared, explained, mu] = pca(features_T);
reduced_data = features_T * coeff(:, 1:2);
figure;
scatter(reduced_data(:, 1), reduced_data(:, 2), 'filled');
text(reduced_data(:, 1)+0.01, reduced_data(:, 2), string(feature_names), 'FontSize', 16);
xlabel('First Principal Component');
ylabel('Second Principal Component');
xlim([-80 80]);
ylim([-80 80]);
title('PCA Projection of Features onto 2D Plane');
grid on;
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\PCA_features.png'))
%% Clusters Wells
runInsiders = false;
if runInsiders
    folderCurrent = folderCurrent + 'Insiders\';
    mkdir(folderCurrent);
end

temp_CAPPatient = string(DataFileWell.CAPPatient(combinedIndex));

cMap = [linspace(0,1,128)'*[1,1], ones(128,1)]; % The blue range
cMap = [cMap; rot90(cMap,2)];                   % add the red range
cMap = (cMap);
    
rowLabels = string(temp_cellType);
for i=1:size(rowLabels, 1)
    rowLabels(i) = rowLabels(i) + '_' + string(i) + '_' + temp_CAPPatient(i);
end
cg = clustergram(normalized_features, ...
                             'ColumnLabels', feature_names,...
                             'RowLabels', rowLabels,...
                             'RowPdist', 'cityblock',...
                             'ColumnPDist', 'cityblock',...
                             'Linkage', 'ward',...
                             'Cluster', 'COLUMN',...
                             'Dendrogram', 0,'Colormap', cMap,...
                           "Standardize", 'column');
%
% GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12; gi13; gi14; gi15; gi16];
% GI = [gi1; gi2; gi5; gi6; gi7; gi8; gi9];
GINums = [6004 5999 5978 5993 5997 5990 5996];
% GINums = [6004 5999 5998 6001 5978];
% GINums = [6004 5975 5976 5989 5954 5980 5984 5992 5978 6001]; % 10
% GINums = [6004 5975 5976 5989 5954 5980 5984 5992 5993 5983 5974 5979 5978]; %good

% GINums = [6004 5975 5976 5964 5960 5954 5980 5984 5952 5971 5985 5987 5983 5974 5979 5978];

%4
% GI = [gi1; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12; gi13; gi14; gi15];
% GINums = [4506 4475 4456 4472 4436 4477 4494 4483 4488 4491 4489 4479 4460 4459];
%3
% GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12; gi13; gi14; gi15; gi16; gi17];
% GINums = [4590 4579 4571 4575 4562 4576 4572 4581 4573 4550 4559 4557 4564 4536 4548 4558 4552];
% GINums = [4590 4579 4571 4575 4562 4576 4572 4581 4573 4550 4559];
%2
% GI = [ng1; ng2; ng3; ng4; ng5; ng6; ng7; ng8; ng9; ng10; gi11; gi12; gi13];
% GINums = [4786 4756 4752 4769 4762 4766 4761 4768 4758 4767 4772 4775 4765];
%1
% GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12];
% GINums = [5927 5909 5914 5911 5903 5912 5898 5905 5907 5913 5910 5900];
% GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8];
% GINums = [6484 6455 6470 6475 6477 6474 6476 6473];
%
% cluster Group
de = 1;
seperate = 0; %33 32 42 41 40

clear groupsCell
numGroups = length(cg.RowLabels) - 1;
numGroups = length(GI);
clusterDE = {};
clusterHC = {};
clusterR = {};
clusterM = {};
clusterS = {};
clusterP = {};
clusterLabels = {};
clusterNames = {};
myCellHC = {};
myCellR = {};
myCellM = {};
myCellS = {};
myCellP = {};
myCellHCCAPS = {};
myCellRCAPS = {};
myCellMCAPS = {};
myCellSCAPS = {};
myCellPCAPS = {};

% groupsCell = cell([numGroups 1]);
groupIndexed=0;
uniqueCAPS = natsort(unique(regexp(string(data.CAPPatient2(combinedIndex)),'\d*','Match', 'once')));
saved_var = cell(1, numGroups);

for groupI=numGroups:-1:1
    group = GINums(groupI);
    myBreakFlag = false;
    flag = false;
    currentGroupInfo = getGroupInfo(cg,group, 1);
    currentGroupInfo = GI(groupI);
    groupsCell(group) = currentGroupInfo;
    if ~runInsiders
        for iterGroup=1:length(clusterLabels)
            currNames = groupsCell(clusterLabels{iterGroup}).RowNodeNames;
            if sum(ismember(currNames, currentGroupInfo.RowNodeNames)) ~= 0
                myBreakFlag = true;
                break;
            end
        end
        if myBreakFlag
           continue; 
        end
    end
    countSevere = sum(contains(string(currentGroupInfo.RowNodeNames), "Severe"));
    countMild = sum(contains(string(currentGroupInfo.RowNodeNames), "Mild"));
    countRisk = sum(contains(string(currentGroupInfo.RowNodeNames), "Premanifest"));
    countHC = sum(contains(string(currentGroupInfo.RowNodeNames), "HC"));
    countProgeria = sum(contains(string(currentGroupInfo.RowNodeNames), "HGPS"));
    
    CAPScores = containers.Map;
    for uC=1:length(uniqueCAPS)
        CAPScores(uniqueCAPS(uC)) = 0;
    end
    inner_values = zeros(1, length(currentGroupInfo.RowNodeNames));
    for item=1:length(currentGroupInfo.RowNodeNames)
        currString = string(currentGroupInfo.RowNodeNames(item));
        pattern = '\d+';
        matches = regexp(currString, pattern, 'match');
        first_number = str2double(matches{1}); % Convert the first match to a number
        inner_values(item) = first_number;

        split_string = strsplit(currString, '_');
        resultCAP = regexp(split_string{3},'\d*','Match', 'once');
        CAPScores(resultCAP) = CAPScores(resultCAP) + 1;
    end
    saved_var{groupI} = inner_values;

    [~,sortedCAPScoresIndexes] = natsort(CAPScores.keys);
    CAPScoresString = string(CAPScores.values);
    CAPScoresString = CAPScoresString(sortedCAPScoresIndexes);
    
    if countMild < countHC && countRisk < countHC && countSevere < countHC && countProgeria < countHC% HC
    if max([countMild countRisk countSevere countProgeria]) / countHC < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countMild countRisk countSevere countProgeria]) / countHC};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    if countMild < countSevere && countRisk < countSevere && countHC < countSevere && countProgeria < countSevere% Severe
    if max([countMild countRisk countHC countProgeria]) / countSevere < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countMild countRisk countHC countProgeria]) / countSevere};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end

    if countSevere < countMild && countRisk < countMild && countHC < countMild && countProgeria < countMild% Mild
    if max([countSevere countRisk countHC countProgeria]) / countMild < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countSevere countRisk countHC countProgeria]) / countMild};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    
    if countSevere < countRisk && countMild < countRisk && countHC < countRisk && countProgeria < countRisk% Risk
    if max([countSevere countMild countHC countProgeria]) / countRisk < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countSevere countMild countHC countProgeria]) / countRisk};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    
    if countSevere < countProgeria && countMild < countProgeria && countHC < countProgeria && countRisk < countProgeria% Progeria
    if max([countSevere countMild countHC countRisk]) / countProgeria < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countSevere countMild countHC countRisk]) / countProgeria};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    
        if flag
            disp("******************************************************");
            checkFlag = false;
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countHC
            checkFlag = true;
           disp(strcat(string(group), "_HC_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "")};
                myCellHC = [myCellHC; strcat('Cluster-',string(groupIndexed), ""), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellHCCAPS = [myCellHCCAPS; strcat('Cluster-',string(groupIndexed), ""), CAPScoresString];

           end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countSevere
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_Severe_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "")};
                myCellS = [myCellS; strcat('Cluster-',string(groupIndexed), ""), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellSCAPS = [myCellSCAPS; strcat('Cluster-',string(groupIndexed), ""), CAPScoresString];
           end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countMild
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_Mild_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_Progeria_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "")};
                myCellM = [myCellM; strcat('Cluster-',string(groupIndexed), ""), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellMCAPS = [myCellMCAPS; strcat('Cluster-',string(groupIndexed), ""), CAPScoresString];

           end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countRisk
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_Risk_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "")};
                myCellR = [myCellR; strcat('Cluster-',string(groupIndexed), ""), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellRCAPS = [myCellRCAPS; strcat('Cluster-',string(groupIndexed), ""), CAPScoresString];
            end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countProgeria
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_HGPS_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "")};
                myCellP = [myCellP; strcat('Cluster-',string(groupIndexed), ""), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellPCAPS = [myCellPCAPS; strcat('Cluster-',string(groupIndexed), ""), CAPScoresString];
            end
        end
        
        disp(strcat("Total_Cells_Is_", string(countMild + countRisk + countSevere + countHC + countProgeria)));

        end
end

clusterNames(1) = {'Cluster 1'};
clusterNames(6) = {'Cluster 2'};
clusterNames(2) = {'Cluster 3'};
clusterNames(3) = {'Cluster 4'};
clusterNames(4) = {'Cluster 5'};
clusterNames(5) = {'Cluster 6'};
clusterNames(7) = {'Cluster 7'};
myCellP(1) = "Cluster 6";
myCellP(2) = "Cluster 7";
myCellR(1) = "Cluster 2";
myCellHC(1) = "Cluster 1";
myCellS(1) = "Cluster 3";
myCellS(2) = "Cluster 4";
myCellS(3) = "Cluster 5";
myCellPCAPS(1) = "Cluster 6";
myCellPCAPS(2) = "Cluster 7";
myCellRCAPS(1) = "Cluster 2";
myCellHCCAPS(1) = "Cluster 1";
myCellSCAPS(1) = "Cluster 3";
myCellSCAPS(2) = "Cluster 4";
myCellSCAPS(3) = "Cluster 5";
groupsCellMotility = groupsCell;
clusterNamesMotility = clusterNames;
%
CM = jet(length(clusterLabels));
CM(:,4) = 0.5;
CMCell = {};
for color=1:size(CM,1)
    CMCell(end+1) = {CM(color,:)};
end
CMCell = CMCell(randperm(numel(CMCell)));
CMCell = {
    [0, 0.2, 0.4],...
    [0.3, 0.8, 1],...
    [1, 0.4, 0.4],...
    [0.8500, 0.3250, 0.0980],...
    [0.5, 0, 0],...
    [0, 0.4470, 0.7410] ,...
    [128, 0, 128]/255,...
};
%
if ~isempty(clusterLabels)
rm = struct('GroupNumber',clusterLabels,'Annotation',clusterNames,...
    'Color',CMCell, 'FontSize', 5);
set(cg,'RowGroupMarker',rm)
cgf = plot(cg);
set(cgf,'FontSize',14, 'fontweight','bold')
fig = get(cgf, 'Parent');
fig.WindowState = 'maximized';
saveas(cgf,strcat(folderCurrent, '\c_Clustergram.png'))
end
%% PCA Features First Component
% features_T = normalized_features;
features_T = z;
% features_T(:,6) = [];
[coeff, score, latent, tsquared, explained] = pca(features_T);
loadings = (coeff(:,1));
[~, most_important_feature_idx] = max(loadings);
most_important_feature = features_T(:, most_important_feature_idx);
disp(most_important_feature_idx);
%%
% Get the loadings for the first principal component
% feature_names_show = feature_names(1:7);
% feature_names_show(6) = [];
loadings = (coeff(:,1));
bar(loadings);
xticks(1:9);
xticklabels(feature_names);
xtickangle(45);
title('Feature Importance');
xlabel('Features');
ylabel('Loadings on First Principal Component');
print('C:\Cell_Migration_MATLAB\IllustratorPlots_4\PCAs\PCA1dbar.png','-dpng','-r900')

%% PCA Features First Two Components
features_T = z;
[coeff, score, latent, tsquared, explained] = pca(features_T);

% Get the absolute value of the loadings for the first two principal components
loadings = abs(coeff(:,1:1));

% Create heatmap
figure;
imagesc(loadings);
colorbar; % shows the color scale
colormap('turbo');
title('Feature Importance');
xlabel('Principal Components');
ylabel('Features');

% Change the x-tick labels to PC1 and PC2
xticks([1]);
xticklabels({'PC1'});

% Change the y-tick labels to feature names
yticks(1:length(feature_names));
yticklabels(feature_names);

% Add text to each cell of the heatmap
for i = 1:size(loadings,1)
    for j = 1:size(loadings,2)
        text(j,i,num2str(loadings(i,j),'%.2f'),...
            'HorizontalAlignment','center',...
            'Color','k');
    end
end
print('C:\Cell_Migration_MATLAB\IllustratorPlots_4\PCAs\PCA1d.png','-dpng','-r900')

%% PCA Features First Two Components
features_T = z;
[coeff, score, latent, tsquared, explained] = pca(features_T);

% Get the absolute value of the loadings for the first two principal components
loadings = abs(coeff(:,1:2));

% Create heatmap
figure;
imagesc(loadings);
colorbar; % shows the color scale
colormap('turbo');
title('Feature Importance');
xlabel('Principal Components');
ylabel('Features');

% Change the x-tick labels to PC1 and PC2
xticks([1 2]);
xticklabels({'PC1', 'PC2'});

% Change the y-tick labels to feature names
yticks(1:length(feature_names));
yticklabels(feature_names);

% Add text to each cell of the heatmap
for i = 1:size(loadings,1)
    for j = 1:size(loadings,2)
        text(j,i,num2str(loadings(i,j),'%.2f'),...
            'HorizontalAlignment','center',...
            'Color','k');
    end
end
print('C:\Cell_Migration_MATLAB\IllustratorPlots_4\PCAs\PCA2d.png','-dpng','-r900')
%%
idxComb = combinedIndex;
xysMain = xys_TRJs_with_HGPS;
xys_full = xys_full_with_HGPS;
for i=1:length(xys_full)
    current = xys_full{i};
    xys{i} = xys_full{i}(:,[3,4]);
end
xys = xys(idxComb);
% xys = xysMain(idxComb);
temp_cellType = DataFileWell.CellType(idxComb);

a1Pol = xys(temp_cellType == "'HC'");
a2Pol = xys(temp_cellType == "'Severe'");
a3Pol = xys(temp_cellType == "'Mild'");
a4Pol = xys(temp_cellType == "'Premanifest'");
a5Pol = xys(temp_cellType == "'HGPS'");

% let's assume classes are stored in a cell array
classes = {a1Pol, a4Pol, a3Pol, a2Pol, a5Pol};
classNames = {'HC', 'Premanifest', 'Mild', 'Severe',  'HGPS'};

figure;
hold on;
maxLag = 20; % Change this value based on your data
colors = ['r', 'g', 'b', 'm', 'c'];

for i = 1:length(classes)
    % get the trajectories for this class
    trajectories = classes{i};
    
    % initialize an array to store the autocorrelation for this class
    autocorrelations = zeros(maxLag+1, length(trajectories));
    
    for j = 1:length(trajectories)
        trajectory = trajectories{j};
        
        % calculate the velocity
        velocity = diff(trajectory);
        
        % calculate the autocorrelation
        [autocorr, lags] = xcorr(velocity, maxLag, 'normalized');
        
        % store the autocorrelation for positive lags only
        autocorrelations(:, j) = autocorr(lags >= 0);
    end
    
    % calculate the mean autocorrelation for this class
    meanAutocorrelation = mean(autocorrelations, 2);
    
    % plot the mean autocorrelation for this class
    plot(0:maxLag, meanAutocorrelation, colors(i));
end

% Set labels and title for the plot
xlabel('Lag');
ylabel('Autocorrelation');
title('Velocity Autocorrelation for Each Class');
legend(classNames);
hold off;


% initialize empty arrays to store all variances and corresponding labels
% allVariances = [];
% allLabels = [];
% 
% for i = 1:length(classes)
%     % get the trajectories for this class
%     trajectories = classes{i};
%     
%     % initialize an array to store the angle variances for this class
%     angleVariances = zeros(length(trajectories), 1);
%     
%     for j = 1:length(trajectories)
%         trajectory = trajectories{j};
%         
%         % calculate the movement vectors
%         movementVectors = diff(trajectory);
%         
%         % calculate the angles between consecutive vectors
%         angles = acos(max(min(dot(movementVectors(1:end-1,:), movementVectors(2:end,:), 2) ./ ...
%             (vecnorm(movementVectors(1:end-1,:), 2, 2) .* vecnorm(movementVectors(2:end,:), 2, 2)), 1), -1));
%         
%         % calculate the variance of the angles
%         angleVariances(j) = var(angles);
%     end
%     
%     % add the variances and labels for this class to the overall arrays
%     allVariances = [allVariances; angleVariances];
%     allLabels = [allLabels; repmat(classNames(i), length(angleVariances), 1)];
% end
% 
% % now, you can create a boxplot for each class
% figure
% boxplot(allVariances, allLabels);
% title('Variance of Turn Angles for each class');


% % Define the bin edges for the histogram (in degrees)
% binEdges = 0:10:180; % Change the bin size according to your requirement
% 
% % Define colors for each class for better differentiation in plot
% colors = ['r', 'g', 'b', 'm', 'c'];
% 
% figure;
% hold on;
% for i = 1:length(classes)
%     % get the trajectories for this class
%     trajectories = classes{i};
%     
%     % initialize an array to store the angle changes for this class
%     angleChanges = [];
%     
%     for j = 1:length(trajectories)
%         trajectory = trajectories{j};
%         
%         % calculate the movement vectors
%         movementVectors = diff(trajectory);
%         
%         % calculate the angles between consecutive vectors
%         angles = acos(max(min(dot(movementVectors(1:end-1,:), movementVectors(2:end,:), 2) ./ ...
%             (vecnorm(movementVectors(1:end-1,:), 2, 2) .* vecnorm(movementVectors(2:end,:), 2, 2)), 1), -1));
%         
%         % append the angles to the angleChanges array
%         angleChanges = [angleChanges; angles];
%     end
%     
%     % convert the angles to degrees
%     angleChanges = rad2deg(angleChanges);
%     
%     % plot the histogram for this class
%     histogram(angleChanges, 'BinEdges', binEdges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', colors(i));
% end
% 
% % Set labels and title for the plot
% xlabel('Angle Change (degrees)');
% ylabel('Probability');
% title('Angular Histogram for Each Class');
% legend(classNames);
% hold off;


% initialize empty arrays to store all aspect ratios and corresponding labels
% allAspectRatios = [];
% allLabels = [];
% 
% for i = 1:length(classes)
%     % get the trajectories for this class
%     trajectories = classes{i};
%     
%     % initialize an array to store the aspect ratio for this class
%     aspectRatios = zeros(length(trajectories), 1);
%     
%     for j = 1:length(trajectories)
%         trajectory = trajectories{j};
%         
%         % compute the aspect ratio
%         aspectRatios(j) = range(trajectory(:, 1)) / range(trajectory(:, 2));
%     end
%     
%     % add the aspect ratios and labels for this class to the overall arrays
%     allAspectRatios = [allAspectRatios; aspectRatios];
%     allLabels = [allLabels; repmat(classNames(i), length(aspectRatios), 1)];
% end
% 
% % now, you can create a boxplot for each class
% figure
% boxplot(allAspectRatios, allLabels);
% title('Aspect ratio of bounding box for each class');


% initialize empty arrays to store all radius of gyrations and corresponding labels
% allRadiusOfGyration = [];
% allLabels = [];
% 
% for i = 1:length(classes)
%     % get the trajectories for this class
%     trajectories = classes{i};
%     
%     % initialize an array to store the radius of gyration for this class
%     radiusOfGyration = zeros(length(trajectories), 1);
%     
%     for j = 1:length(trajectories)
%         trajectory = trajectories{j};
%         x = trajectory(:, 1) - mean(trajectory(:, 1));
%         y = trajectory(:, 2) - mean(trajectory(:, 2));
%         
%         % compute the radius of gyration
%         radiusOfGyration(j) = sqrt(mean(x.^2 + y.^2));
%     end
%     
%     % add the radius of gyration and labels for this class to the overall arrays
%     allRadiusOfGyration = [allRadiusOfGyration; radiusOfGyration];
%     allLabels = [allLabels; repmat(classNames(i), length(radiusOfGyration), 1)];
% end
% 
% % now, you can create a boxplot for each class
% figure
% boxplot(allRadiusOfGyration, allLabels);
% title('Radius of gyration for each class');

% initialize an empty cell array to store dominant frequencies for each class
% dominantFrequencies = cell(1, length(classes));
% 
% for i = 1:length(classes)
%     % get the trajectories for this class
%     trajectories = classes{i};
%     
%     % initialize an array to store the dominant frequencies for this class
%     dominantFrequencies{i} = zeros(length(trajectories), 1);
%     
%     for j = 1:length(trajectories)
%         trajectory = trajectories{j};
%         x = trajectory(:, 1) - mean(trajectory(:, 1));
%         y = trajectory(:, 2) - mean(trajectory(:, 2));
%         angles = atan2(y, x);
%         Y = fft(angles);
%         Pyy = Y.*conj(Y)/length(angles);
%         
%         % compute the dominant frequency
%         [~, index] = max(Pyy);
%         dominantFrequencies{i}(j) = index;
%     end
% end
% 
% % now, you can create a histogram for each class
% for i = 1:length(classes)
%     subplot(length(classes), 1, i);
%     histogram(dominantFrequencies{i});
%     title(classNames{i});
%     xlim([0 10]); 
%     ylim([0 1000]); 
% end
%% Clusters dTheta and angular velocity
saved_var_flipped = flip(saved_var);
num_clusters = length(saved_var_flipped);
xys = cell(length(xys_full),1);
for i=1:length(xys_full)
    current = xys_full{i};
    xys{i} = xys_full{i}(:,[3,4]);
end
xys = xys(idxComb);
% xys = xysMain(idxComb);

clusters_order = [1,6,2,3,4,5,7];
clusters_order1 = [1,3,4,5,6,2,7];

for i = 1:num_clusters
    comb = saved_var_flipped{i}';
    aPol = xys(comb);
    % dTheta PDF
    close all
    param.showfig=1;
    param.saveres=0;
    param.dim=2;
    param.outfigurenum=405;
    param.binnum=6;
    param.markertype='-';
    tloi = [1 2 3 4 5 6];
    cid=copper(length(tloi));
    % param.markertype=newcolorsRGB{1};
    % figure();
    get_dtheta_PDF(aPol, tloi, param);
    % legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
    title(strcat('Cluster','{ }', string(clusters_order1(i)) ,' dθ PDF'));
    xlabel('dθ') ;
    ylabel('Occurrence');
        axis square;
        box on
            colorMatrix = repmat(linspace(1, length(tloi), 100), length(tloi), 1);
            color_ax = axes('Position', [0.7, 0.1, 0.2, 0.05]);
            contourf(color_ax, colorMatrix, length(tloi), 'LineStyle', 'none');
            colormap(color_ax, cid);
            axis off;
            color_bar_title = title(color_ax, '0 - 60 Minutes', 'Units', 'normalized', 'Position', [0.5, 1.3, 0],'FontSize', 18, 'FontWeight', 'bold');

    g = gcf;
    g.WindowState = 'maximized';
    saveas(g, strcat(folderCurrent, '\clusters_z_Theta_PDF_',string(clusters_order1(i)),'.png'))
    
    
end
%
    
param.showfig=1;
param.saveres=0;
param.dim=2;
param.outfigurenum=305;
param.binnum=35;
CMCell = {
    [0, 0.2, 0.4],...
    [0.3, 0.8, 1],...
    [1, 0.4, 0.4],...
    [0.8500, 0.3250, 0.0980],...
    [0.5, 0, 0],...
    [0, 0.4470, 0.7410] ,...
    [128, 0, 128]/255,...
};
for i = 1:num_clusters
    comb = saved_var_flipped{i}';
    aPol = xys(comb);

    param.markertype=CMCell{i};
    get_dR_polarity(aPol, 1, param);
    hold on;
    
    
end
title('Clusters dR Polarity');

g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\clusters_z_Theta.png'))

param.showfig=1;
param.saveres=0;
param.markertype='.-';
param.outfigurenum=301;
param.dim=2;
param.linear=0;
param.aviv=1;

for j = 1:num_clusters
    i = clusters_order(j);
    comb = saved_var_flipped{i}';
    aPol = xys(comb);

    param.markertype=CMCell{i};
    get_MSD(aPol, 10, param);
    hold on;
    
    
end

param.markertype='k:';
param.linear=1;
get_MSD(aPol, 10, param);
hold off;
title('Mean Squared Displacement');

ylabel('MSD (µm2)') ;
xlabel('Time lag (min)');

legend('Cluster 1','Cluster 2' ,'Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7', 'Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\clusters_z_MSD.png'))
%%
saved_var_flipped = flip(saved_var);
num_clusters = length(saved_var_flipped);
clusters_order = [1,6,2,3,4,5,7];
clusters_order1 = [1,3,4,5,6,2,7];

xys = cell(length(xys_full),1);
for i=1:length(xys_full)
    current = xys_full{i};
    xys{i} = xys_full{i}(:,[3,4]);
end
xys = xys(idxComb);

param.dim=2; 
param.tlag=1;
param.saveres=0;
param.showfig=1;       
param.markertype2='none'; % HC
param.outfigurenum=303;   
param.aviv=1;
param.aviv2 = 0;
param.aviv3 = 1;

for j = 1:num_clusters
        i = clusters_order(j);
    comb = saved_var_flipped{i}';
    aPol = xys(comb);

    param.markertype=CMCell{i};
    get_ACF1(aPol, 10, param);
    hold on;
    
    
end

title('Auto-Correlation Function');
ylabel('ACF') ;
xlabel('Time lag (min)');
% legend('HC','line','Premanifest' ,'line','Mild', 'line','Severe','line','HGPS','line','Location','bestoutside');
legend('Cluster 1','Cluster 2' ,'Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7', 'Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\clusters_z_ACF.png'))


%% dR PDF
% close all force
param.dxmax=50;
param.binn=70;            
param.showfig=1;
param.saveres=0;
param.dim=2;            
param.outfigurenum=302;
param.aviv=0;

param.markertype2=':';
for j = 1:num_clusters
        i = clusters_order(j);
    comb = saved_var_flipped{i}';
    aPol = xys(comb);

    param.markertype=CMCell{i};
    get_dR_PDF(aPol, 10, param);
    hold on;
    
    
end

title('PDF Cellular Displacements');
ylabel('Occurrence') ;
xlabel('Displacement (µm)');
legend('Cluster 1','Cluster 2' ,'Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7', 'Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\clusters_z_PDF.png'))
%% Draw (project) 1d of the 2dPCA clusters line
temp_Patient;
new_feat = round(temp_Patient);
patients;
% array([  0.,  60.,  69.,  75.,  76.,  76.,  84.,  91.,  94.,  96., 100.,
%        104., 106., 110. 13, 111., 111., 114., 115.17, 126., 126., 129., 131.,
%        132., 141., 999.]) 0.992910311002455 is interpolation for matlab
%        119
map = [1, 0.9628962515657298, 0.8832500485011887, 0.9315501070412957, 1.0475624836070825, ...
    0.9062262706482912, 0.9386139949131951, 0.9876259986191774, 1.001472534625067, 0.9328344061142667, ...
    0.9464954566015233, 0.9367413159630704, 0.9382691582210567, 0.9748999561999504, 1.002165056468319, ...
    0.9692848270579759, 0.992910311002455, 1.016535794946934, 0.9373562717737963, 1.0408265727672568 ...
    1.272197713272603]';
uP = round(unique(temp_Patient));
for tP=1:length(uP)
    new_feat(round(temp_Patient) == uP(tP)) = map(tP);
end
features_T = z;
[coeff, score, latent, tsquared, explained, mu] = pca(features_T);
reduced_data = features_T * coeff(:, 1:2);
saved_var_flipped = flip(saved_var);
num_clusters = length(saved_var_flipped);
num_points = length(reduced_data);
IDX = zeros(num_points, 1);

p = polyfit(reduced_data(:,1), reduced_data(:,2), 1);
m = p(1);
b = p(2);
x = reduced_data(:,1);
y = reduced_data(:,2);
x_projected = (x + m*y - m*b) ./ (m^2 + 1);

% Assign cluster numbers to the data points
for i = 1:num_clusters
    IDX(saved_var_flipped{i}) = i;
end

% Create the PCA plot
figure;
hold on;

for i = 1:num_clusters - 1
    cluster_points = x_projected(IDX == i);
    cluster_patients = temp_Patient(IDX == i);
%     cluster_features = features_T(IDX == i, 4);
    cluster_features = new_feat(IDX == i);

    
    % Scatter plot for each cluster
    scatter(cluster_points, cluster_features, 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
  
end
for i = 1:num_clusters - 1
    cluster_points = x_projected(IDX == i);
    cluster_patients = temp_Patient(IDX == i);
%     cluster_features = features_T(IDX == i, 4);
    cluster_features = new_feat(IDX == i);
    
    % Compute the mean for each cluster
    mean_cluster_point = mean(cluster_points);
    mean_cluster_feature = mean(cluster_features);
    
    % Offset to position the star and text above the data points
    y_offset = 0.05 * (max(cluster_features) - min(cluster_features)); % Adjust as needed
    
    % Plot the black star at the mean
    scatter(mean_cluster_point, mean_cluster_feature + y_offset, 100, 'k', '*'); % 100 is the marker size; adjust as necessary
    
    % Annotate the star with the cluster number
    text(mean_cluster_point, mean_cluster_feature + y_offset, sprintf(clusterNames{i}), ...
        'Color', 'k', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
end

% Adding legend to the plot
legend(string(clusterNames(1:end-1)), 'Location', 'best'); % or 'Location', 'northeast', or wherever best fits your data
xlabel('Projected PCA');
ylabel('Feature Value');
title('Scatter plot of projected PCA vs. feature');
grid on;

hold off;
%% PCA Groups
temp_Patient;
new_feat = round(temp_Patient);
patients;
% array([  0.,  60.,  69.,  75.,  76.,  76.,  84.,  91.,  94.,  96., 100.,
%        104., 106., 110. 13, 111., 111., 114., 115.17, 126., 126., 129., 131.,
%        132., 141., 999.]) 0.992910311002455 is interpolation for matlab
%        119
% map = [1, 0.9628962515657298, 0.8832500485011887, 0.9315501070412957, 1.0475624836070825, ...
%     0.9062262706482912, 0.9386139949131951, 0.9876259986191774, 1.001472534625067, 0.9328344061142667, ...
%     0.9464954566015233, 0.9367413159630704, 0.9382691582210567, 0.9748999561999504, 1.002165056468319, ...
%     0.9692848270579759, 0.992910311002455, 1.016535794946934, 0.9373562717737963, 1.0408265727672568 ...
%     1.272197713272603]'; % widthToLength
map = [1, 1.0806452120803148, 1.0739396305528609, 1.083671045294664, 0.9312406164134382, ...
    0.8751886409111126, 1.034898150435135, 1.1001504787102776, 1.1578667759821568, 1.160018265675352, ...
    0.8284036633274416, 1.1589048578799673, 0.7034004682316395, 0.70351551967129945, 1.111981913103297, ...
    1.2910475879181663, 1.4516386057056714, 1.5575123647002798, 1.6633861236948882, 1.7772495419577079, ...
    0.3475618917471596]'; % LaminB1
% mapHGPS = 0.3475618917471596; % LaminB1
uP = round(unique(temp_Patient));
for tP=1:length(uP)
    new_feat(round(temp_Patient) == uP(tP)) = map(tP);
end
% new_feat(new_feat < 0.4) = 1;
% Create the colormap from light green to dark green
num_colors = length(map) - 1; % Or some other appropriate number for granularity
colormap_green = [linspace(0.8, 0, num_colors)' linspace(1, 0.5, num_colors)' linspace(0.8, 0, num_colors)'];
colormap_purple = [linspace(0.8, 0.4, num_colors)' linspace(0.8, 0, num_colors)' linspace(1, 0.5, num_colors)'];
uF = unique(new_feat);
min_val = min(uF(1:end));
max_val = max(uF(1:end));
new_feat_normalized = (new_feat - min_val) / (max_val - min_val);
% Map normalized new_feat values to the colormap
color_indices = ceil(new_feat_normalized * (num_colors-1)) + 1;

colormap_use = colormap_purple;

colors_mapped = colormap_use(color_indices, :);


features_T = z;
[coeff, score, latent, tsquared, explained, mu] = pca(features_T);
reduced_data = features_T * coeff(:, 1:2);
num_clusters = 5;
num_points = length(reduced_data);
IDX = zeros(num_points, 1);
markers = {'^', 'o', '>', 'd', 's', 'h', 'v', '<'};

% Assign cluster numbers to the data points
classNames = {'HC', 'Premanifest', 'Mild', 'Severe',  'HGPS'};

temp_cellType = DataFileWell.CellType(idxComb);
IDX(temp_cellType == "'HC'") = 1;
IDX(temp_cellType == "'Severe'") = 4;
IDX(temp_cellType == "'Mild'") = 3;
IDX(temp_cellType == "'Premanifest'") = 2;
IDX(temp_cellType == "'HGPS'") = 5;

% Create the PCA plot
figure;
hold on;

colors = lines(num_clusters); % Creates a color matrix for the clusters
legend_entries = cell(1, num_clusters);
for i = num_clusters:-1:1
    cluster_points = reduced_data(IDX == i, :);
    % Get the corresponding colors for this cluster
    if i == num_clusters
        cluster_colors = colors_mapped(IDX == i, :);
            scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
    else
    cluster_colors = colors_mapped(IDX == i, :);
    if i == 1 || i == 2
        o = 0.15;
    else
        o = 0.75;
    end
        scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha',o, 'MarkerEdgeAlpha', o);
    end

    
end
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    mean_point = mean(cluster_points, 1);
    co = 'k';
    if i == num_clusters
        co = 'k';
    end
    s2 = scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerEdgeColor', co, 'LineWidth', 1.7, 'Marker', 'p', 'SizeData', 100);
    text(mean_point(1)+0.05, mean_point(2), string(classNames(i)), 'FontSize', 12);
end
xlabel('First Principal Component');
ylabel('Second Principal Component');
xlim([-35 5]);
ylim([-15 20]);
title('PCA Projection of Groups onto 2D Plane');
colorbar;
colormap(colormap_use);
caxis([min_val, max_val]);
% Add text next to the colorbar:
ax = gca;
% text(ax.Position(1) + ax.Position(3) + 8, ax.Position(2) + ax.Position(4) / 2, ...
%      'Cells Width to Length Ratio', 'Rotation', 90, 'VerticalAlignment', 'middle');
text(ax.Position(1) + ax.Position(3) + 8, ax.Position(2) + ax.Position(4) / 2, ...
     'LaminB1 Total Area', 'Rotation', 90, 'VerticalAlignment', 'middle');
grid on;
hold off;

legend(string(classNames), 'Location', 'northeast');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
% exportgraphics(g, strcat(folderCurrent, '\PCA2d_Groups_ratio.png'), 'Resolution', 900)
exportgraphics(g, strcat(folderCurrent, '\PCA2d_Groups_LaminB1.png'), 'Resolution', 900)


figure;
hold on;

colors = lines(num_clusters); % Creates a color matrix for the clusters
legend_entries = cell(1, num_clusters);
for i = num_clusters:-1:1
    cluster_points = reduced_data(IDX == i, :);
    if i == num_clusters
        cluster_colors = colors_mapped(IDX == i, :);
    scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
    else
    cluster_colors = colors_mapped(IDX == i, :);
    if i == 1 || i == 2
        o = 0.35;
    else
        o = 0.75;
    end
        scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha',o, 'MarkerEdgeAlpha', o);
    end



end
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    mean_point = mean(cluster_points, 1);
    co = 'k';
    if i == num_clusters
        co = 'k';
    end
    s2 = scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerEdgeColor', co, 'LineWidth', 2.5, 'Marker', 'p', 'SizeData', 100);
        text(mean_point(1)+0.05, mean_point(2), string(classNames(i)), 'FontSize', 16,'fontweight','bold');
end

xlabel('First Principal Component');
ylabel('Second Principal Component');
xlim([-6 2.5]);
ylim([-2 4]);
title('Close-up of the PCA Projection of Groups onto 2D Plane');
colorbar;
colormap(colormap_use);
caxis([min_val, max_val]);
ax = gca;
% text(ax.Position(1) + ax.Position(3) + 2.5, ax.Position(2) + ax.Position(4) / 2, ...
%      'Cells Width to Length Ratio', 'Rotation', 90, 'VerticalAlignment', 'middle');
text(ax.Position(1) + ax.Position(3) + 2.5, ax.Position(2) + ax.Position(4) / 2, ...
     'LaminB1 Total Area', 'Rotation', 90, 'VerticalAlignment', 'middle');
grid on;
hold off;


legend(string(classNames), 'Location', 'best');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
% exportgraphics(g, strcat(folderCurrent, '\PCA2d_Groups_Closeup_ratio.png'), 'Resolution', 900)
exportgraphics(g, strcat(folderCurrent, '\PCA2d_Groups_Closeup_LaminB1.png'), 'Resolution', 900)

%% PCA Clusters

clusterNames(1) = {'Cluster 1'};
clusterNames(6) = {'Cluster 2'};
clusterNames(2) = {'Cluster 3'};
clusterNames(3) = {'Cluster 4'};
clusterNames(4) = {'Cluster 5'};
clusterNames(5) = {'Cluster 6'};
clusterNames(7) = {'Cluster 7'};

temp_Patient;
new_feat = round(temp_Patient);
patients;
% array([  0.,  60.,  69.,  75.,  76.,  76.,  84.,  91.,  94.,  96., 100.,
%        104., 106., 110. 13, 111., 111., 114., 115.17, 126., 126., 129., 131.,
%        132., 141., 999.]) 0.992910311002455 is interpolation for matlab
%        119
map = [1, 0.9628962515657298, 0.8832500485011887, 0.9315501070412957, 1.0475624836070825, ...
    0.9062262706482912, 0.9386139949131951, 0.9876259986191774, 1.001472534625067, 0.9328344061142667, ...
    0.9464954566015233, 0.9367413159630704, 0.9382691582210567, 0.9748999561999504, 1.002165056468319, ...
    0.9692848270579759, 0.992910311002455, 1.016535794946934, 0.9373562717737963, 1.0408265727672568 ...
    1.272197713272603]'; % widthToLength
% map = [1, 1.0806452120803148, 1.0739396305528609, 1.083671045294664, 0.9312406164134382, ...
%     0.8751886409111126, 1.034898150435135, 1.1001504787102776, 1.1578667759821568, 1.160018265675352, ...
%     0.8284036633274416, 1.1589048578799673, 0.7034004682316395, 0.70351551967129945, 1.111981913103297, ...
%     1.2910475879181663, 1.4516386057056714, 1.5575123647002798, 1.6633861236948882, 1.7772495419577079, ...
%     0.3475618917471596]'; % LaminB1
% mapHGPS = 0.3475618917471596; % LaminB1
uP = round(unique(temp_Patient));
for tP=1:length(uP)
    new_feat(round(temp_Patient) == uP(tP)) = map(tP);
end
% new_feat(new_feat < 0.4) = 1;
% Create the colormap from light green to dark green
num_colors = length(map) - 1; % Or some other appropriate number for granularity
colormap_green = [linspace(0.8, 0, num_colors)' linspace(1, 0.5, num_colors)' linspace(0.8, 0, num_colors)'];
colormap_purple = [linspace(0.8, 0.4, num_colors)' linspace(0.8, 0, num_colors)' linspace(1, 0.5, num_colors)'];
colormap_red = [linspace(1, 0.8, num_colors)' linspace(0.8, 0, num_colors)' linspace(0.8, 0, num_colors)'];
colormap_blue = [linspace(0.8, 0, num_colors)' linspace(0.8, 0, num_colors)' linspace(1, 0.8, num_colors)'];
uF = unique(new_feat);
min_val = min(uF(1:end));
max_val = max(uF(1:end));
new_feat_normalized = (new_feat - min_val) / (max_val - min_val);
% Map normalized new_feat values to the colormap
color_indices = ceil(new_feat_normalized * (num_colors-1)) + 1;

colormap_use = colormap_blue;

colors_mapped = colormap_use(color_indices, :);


features_T = z;
[coeff, score, latent, tsquared, explained, mu] = pca(features_T);
reduced_data = features_T * coeff(:, 1:2);
saved_var_flipped = flip(saved_var);
num_clusters = length(saved_var_flipped);
num_points = length(reduced_data);
IDX = zeros(num_points, 1);
markers = {'^', 'o', '>', 'd', 's', 'h', 'v', '<'};

% Assign cluster numbers to the data points
for i = 1:num_clusters
    IDX(saved_var_flipped{i}) = i;
end

% Create the PCA plot
figure;
hold on;

colors = lines(num_clusters); % Creates a color matrix for the clusters
legend_entries = cell(1, num_clusters);
for i = num_clusters:-1:1
    cluster_points = reduced_data(IDX == i, :);
    % Get the corresponding colors for this cluster
    if i == num_clusters
        cluster_colors = colors_mapped(IDX == i, :);
            scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
    else
    if i == 1 || i == 2
        o = 0.75;
    else
        o = 0.75;
    end
    cluster_colors = colors_mapped(IDX == i, :);
        scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha', o, 'MarkerEdgeAlpha', o);
    end
    
end
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    mean_point = median(cluster_points, 1);
    co = 'white';
    if i == num_clusters
        co = 'k';
    end
    s2 = scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerEdgeColor', co, 'LineWidth', 1.7, 'Marker', 'p', 'SizeData', 100);
    text(mean_point(1)+0.05, mean_point(2), string(clusterNames(i)), 'FontSize', 12);
end

for i = 1:num_clusters
    representative_color = colors_mapped(IDX == i, :);
    representative_color = representative_color(1, :);
    p(i) = scatter(NaN, NaN, 25, representative_color , markers{i}, 'filled');
end

xlabel('First Principal Component');
ylabel('Second Principal Component');
xlim([-35 5]);
ylim([-15 20]);
title('PCA Projection of Clusters onto 2D Plane');
colorbar;
colormap(colormap_use);
caxis([min_val, max_val]);
% Add text next to the colorbar:
ax = gca;
text(ax.Position(1) + ax.Position(3) + 8, ax.Position(2) + ax.Position(4) / 2, ...
     'Cells Width to Length Ratio', 'Rotation', 90, 'VerticalAlignment', 'middle');
% text(ax.Position(1) + ax.Position(3) + 8, ax.Position(2) + ax.Position(4) / 2, ...
%      'LaminB1 Total Area', 'Rotation', 90, 'VerticalAlignment', 'middle');
grid on;
hold off;

% Add a legend
newClusterNames = cell(1, numel(clusterNames) * 2);
newString = 'NewString';
for i = 1:numel(clusterNames)
    newClusterNames{2 * i - 1} = clusterNames{i};
    
    if i < numel(clusterNames)
        newClusterNames{2 * i} = strcat("Center of ",clusterNames{i});
    end
end
newClusterNames{end} = strcat("Center of ",clusterNames{end});
newOrder = [1 6 2 3 4 5 7];
p = p(newOrder);
clusterNames = clusterNames(newOrder);

lgd = legend(p, string(clusterNames), 'Location', 'bestoutside');
lgd.FontSize = 16;
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
exportgraphics(g, strcat(folderCurrent, '\PCA2d_Clusters_ratio_g.png'), 'Resolution', 900)
% exportgraphics(g, strcat(folderCurrent, '\PCA2d_Clusters_LaminB1.png'), 'Resolution', 900)

clusterNames(1) = {'Cluster 1'};
clusterNames(6) = {'Cluster 2'};
clusterNames(2) = {'Cluster 3'};
clusterNames(3) = {'Cluster 4'};
clusterNames(4) = {'Cluster 5'};
clusterNames(5) = {'Cluster 6'};
clusterNames(7) = {'Cluster 7'};
figure;
hold on;

colors = lines(num_clusters); % Creates a color matrix for the clusters
legend_entries = cell(1, num_clusters);
for i = num_clusters:-1:1
    cluster_points = reduced_data(IDX == i, :);
    if i == num_clusters
        cluster_colors = colors_mapped(IDX == i, :);
            scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
    else
    if i == 1 || i == 2
        o = 0.35;
    else
        o = 0.75;
    end
    cluster_colors = colors_mapped(IDX == i, :);
        scatter(cluster_points(:, 1), cluster_points(:, 2), 25, cluster_colors, markers{i}, 'filled', ...
        'MarkerFaceAlpha', o, 'MarkerEdgeAlpha', o);
    end


end
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    mean_point = median(cluster_points, 1);
    co = 'k';
    if i == num_clusters
        co = 'k';
    end
    s2 = scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerEdgeColor', co, 'LineWidth', 2.5, 'Marker', 'p', 'SizeData', 100);
        text(mean_point(1)+0.05, mean_point(2), string(clusterNames(i)), 'FontSize', 30,'fontweight','bold');
end
for i = 1:num_clusters
    representative_color = colors_mapped(IDX == i, :);
    representative_color = representative_color(1, :);
    p(i) = scatter(NaN, NaN, 25, representative_color , markers{i}, 'filled');
end
xlabel('First Principal Component');
ylabel('Second Principal Component');
xlim([-6 3.5]);
ylim([-2 4]);
title('Close-up of the PCA Projection of Clusters onto 2D Plane');
colorbar;
colormap(colormap_use);
caxis([min_val, max_val]);
% ax = gca;
text(ax.Position(1) + ax.Position(3) + 3.5, ax.Position(2) + ax.Position(4) / 2, ...
     'Cells Width to Length Ratio', 'Rotation', 90, 'VerticalAlignment', 'middle');
%  text(ax.Position(1) + ax.Position(3) + 3.5, ax.Position(2) + ax.Position(4) / 2, ...
%      'LaminB1 Total Area', 'Rotation', 90, 'VerticalAlignment', 'middle');
grid on;
hold off;

% Add a legend
newClusterNames = cell(1, numel(clusterNames) * 2);
newString = 'NewString';
for i = 1:numel(clusterNames)
    newClusterNames{2 * i - 1} = clusterNames{i};
    
    if i < numel(clusterNames)
        newClusterNames{2 * i} = strcat("Center of ",clusterNames{i});
    end
end
newClusterNames{end} = strcat("Center of ",clusterNames{end});
newOrder = [1 6 2 3 4 5 7];
p = p(newOrder);
clusterNames = clusterNames(newOrder);

lgd = legend(p, string(clusterNames), 'Location', 'bestoutside');
lgd.FontSize = 16;
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
exportgraphics(g, strcat(folderCurrent, '\PCA2d_Clsuters_Closeup_ratio_g.png'), 'Resolution', 900)
% exportgraphics(g, strcat(folderCurrent, '\PCA2d_Clsuters_Closeup_LaminB1.png'), 'Resolution', 900)

%% PCA Clsuters
clusterNames(1) = {'Cluster 1'};
clusterNames(6) = {'Cluster 2'};
clusterNames(2) = {'Cluster 3'};
clusterNames(3) = {'Cluster 4'};
clusterNames(4) = {'Cluster 5'};
clusterNames(5) = {'Cluster 6'};
clusterNames(7) = {'Cluster 7'};

features_T = z;
[coeff, score, latent, tsquared, explained, mu] = pca(features_T);
reduced_data = features_T * coeff(:, 1:2);
saved_var_flipped = flip(saved_var);
num_clusters = length(saved_var_flipped);
num_points = length(reduced_data);
IDX = zeros(num_points, 1);

% Assign cluster numbers to the data points
for i = 1:num_clusters
    IDX(saved_var_flipped{i}) = i;
end

% Create the PCA plot
figure;
hold on;

colors = lines(num_clusters); % Creates a color matrix for the clusters
legend_entries = cell(1, num_clusters);
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    s1 = scatter(cluster_points(:, 1), cluster_points(:, 2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    mean_point = median(cluster_points, 1);
    uistack(s1, 'bottom');
    co = 'white';
    if i == num_clusters
        co = 'k';
    end
    s2 = scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerEdgeColor', co, 'LineWidth', 1.7, 'Marker', 'p', 'SizeData', 100);
    text(mean_point(1)+0.05, mean_point(2), string(clusterNames(i)), 'FontSize', 12);
    uistack(s2, 'top');
end
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    p(i) = scatter(NaN, NaN, 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)));
end
xlabel('First Principal Component');
ylabel('Second Principal Component');
xlim([-35 5]);
ylim([-15 20]);
title('PCA Projection of Clusters onto 2D Plane');
grid on;
hold off;

% Add a legend
newClusterNames = cell(1, numel(clusterNames) * 2);
newString = 'NewString';
for i = 1:numel(clusterNames)
    newClusterNames{2 * i - 1} = clusterNames{i};
    
    if i < numel(clusterNames)
        newClusterNames{2 * i} = strcat("Center of ",clusterNames{i});
    end
end
newClusterNames{end} = strcat("Center of ",clusterNames{end});
newOrder = [1 6 2 3 4 5 7];
p = p(newOrder);
clusterNames = clusterNames(newOrder);

lgd = legend(p,string(clusterNames), 'Location', 'bestoutside');
lgd.FontSize = 16;

    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
exportgraphics(g, strcat(folderCurrent, '\PCA2d.png'), 'Resolution', 900)

clusterNames(1) = {'Cluster 1'};
clusterNames(6) = {'Cluster 2'};
clusterNames(2) = {'Cluster 3'};
clusterNames(3) = {'Cluster 4'};
clusterNames(4) = {'Cluster 5'};
clusterNames(5) = {'Cluster 6'};
clusterNames(7) = {'Cluster 7'};
figure;
hold on;

colors = lines(num_clusters); % Creates a color matrix for the clusters
legend_entries = cell(1, num_clusters);
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    s1 = scatter(cluster_points(:, 1), cluster_points(:, 2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
    mean_point = median(cluster_points, 1);
    uistack(s1, 'bottom');
    co = 'k';
    if i == num_clusters
        co = 'k';
    end
    s2 = scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), ...
        'MarkerEdgeColor', co, 'LineWidth', 2.5, 'Marker', 'p', 'SizeData', 100);
        text(mean_point(1)+0.05, mean_point(2), string(clusterNames(i)), 'FontSize', 30,'fontweight','bold');
    uistack(s2, 'top');
end
for i = 1:num_clusters
    cluster_points = reduced_data(IDX == i, :);
    p(i) = scatter(NaN, NaN, 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)));
end
xlabel('First Principal Component');
ylabel('Second Principal Component');
xlim([-6 3.5]);
ylim([-2 4]);
title('Close-up of the PCA Projection of Clusters onto 2D Plane');
grid on;
hold off;

% Add a legend
newClusterNames = cell(1, numel(clusterNames) * 2);
newString = 'NewString';
for i = 1:numel(clusterNames)
    newClusterNames{2 * i - 1} = clusterNames{i};
    
    if i < numel(clusterNames)
        newClusterNames{2 * i} = strcat("Center of ",clusterNames{i});
    end
end
newClusterNames{end} = strcat("Center of ",clusterNames{end});
newOrder = [1 6 2 3 4 5 7];
p = p(newOrder);
clusterNames = clusterNames(newOrder);

lgd = legend(p, string(clusterNames), 'Location', 'bestoutside');
lgd.FontSize = 16;
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
exportgraphics(g, strcat(folderCurrent, '\PCA2d_Closeup.png'), 'Resolution', 900)
%%
% reduced_data_tsne = tsne(features_T);
% 
% % Create IDX vector
% num_clusters = length(saved_var_flipped);
% num_points = length(reduced_data_tsne);
% IDX = zeros(num_points, 1);
% 
% % Assign cluster numbers to the data points
% for i = 1:num_clusters
%     IDX(saved_var_flipped{i}) = i;
% end
% 
% % Create the t-SNE plot
% % figure;
% % hold on;
% 
% colors = lines(num_clusters); % Creates a color matrix for the clusters
% legend_entries = cell(1, num_clusters);
% 
% % for i = 1:num_clusters
% %     cluster_points = reduced_data_tsne(IDX == i, :);
% %     scatter(cluster_points(:, 1), cluster_points(:, 2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)));
% %     legend_entries{i} = sprintf('Cluster %d', i);
% % end
% % 
% % xlabel('t-SNE Dimension 1');
% % ylabel('t-SNE Dimension 2');
% % title('t-SNE Projection of Clusters');
% % grid on;
% % hold off;
% 
% % Add a legend
% % legend(string(clusterNames), 'Location', 'best');
% % g = gcf;
% % g.WindowState = 'maximized';
% % saveas(g, strcat(folderCurrent, '\TSNE2d.png'))
% % figure;
% % scatter(reduced_data(:, 1), reduced_data(:, 2), 'filled');
% % % text(reduced_data(:, 1)+0.01, reduced_data(:, 2), string(feature_names), 'FontSize', 16);
% % xlabel('First Principal Component');
% % ylabel('Second Principal Component');
% % xlim([-35 5]);
% % ylim([-15 20]);
% % title('PCA Projection of Clusters onto 2D Plane');
% % grid on;
% num_clusters = length(saved_var_flipped);
% num_points = length(reduced_data);
% IDX = zeros(num_points, 1);
% 
% % Assign cluster numbers to the data points
% for i = 1:num_clusters
%     IDX(saved_var_flipped{i}) = i;
% end
% 
% % Create the 2D PCA plot with cluster mean points
% % figure;
% % hold on;
% 
% colors = lines(num_clusters); % Creates a color matrix for the clusters
% legend_entries = cell(1, num_clusters);
% 
% % for i = 1:num_clusters
% %     cluster_points = reduced_data(IDX == i, :);
% %     mean_point = mean(cluster_points, 1);
% %     scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% %     legend_entries{i} = sprintf('Cluster %d', i);
% % end
% % 
% % xlabel('First Principal Component');
% % ylabel('Second Principal Component');
% % title('PCA Projection of Cluster Mean Points onto 2D Plane');
% % grid on;
% % hold off;
% % legend(string(clusterNames), 'Location', 'best');
% % g = gcf;
% % g.WindowState = 'maximized';
% % saveas(g, strcat(folderCurrent, '\PCA2d_2.png'))
% 
% num_clusters = length(saved_var_flipped);
% num_points = length(reduced_data_tsne);
% IDX = zeros(num_points, 1);
% 
% % Assign cluster numbers to the data points
% for i = 1:num_clusters
%     IDX(saved_var_flipped{i}) = i;
% end
% 
% % Create the 2D t-SNE plot with cluster mean points
% % figure;
% % hold on;
% 
% colors = lines(num_clusters); % Creates a color matrix for the clusters
% legend_entries = cell(1, num_clusters);
% 
% % for i = 1:num_clusters
% %     cluster_points = reduced_data_tsne(IDX == i, :);
% %     mean_point = mean(cluster_points, 1);
% %     scatter(mean_point(1), mean_point(2), 'filled', 'MarkerFaceColor', cell2mat(CMCell(i)), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% %     legend_entries{i} = sprintf('Cluster %d', i);
% % end
% % 
% % xlabel('t-SNE Dimension 1');
% % ylabel('t-SNE Dimension 2');
% % title('t-SNE Projection of Cluster Mean Points onto 2D Plane');
% % grid on;
% % hold off;
% % legend(string(clusterNames), 'Location', 'best');
% % g = gcf;
% % g.WindowState = 'maximized';
% % saveas(g, strcat(folderCurrent, '\TSNE2d_2.png'))
%
% close all force

%$%$ BAR PLOT
myCellTotal = [myCellHC; myCellR; myCellM; myCellS; myCellP];
XAxis = categorical(myCellTotal(:,1));
XAxis = reordercats(XAxis,flip(myCellTotal(:,1)));
YAxis = str2double(myCellTotal(:,2:end));
sumY = sum(YAxis, 2);
Answer = bsxfun(@rdivide,YAxis,sumY(:));

YAxis1 = Answer(:,1);
YAxis2 = Answer(:,2);
YAxis3 = Answer(:,3);
YAxis4 = Answer(:,4);
YAxis5 = Answer(:,5);

Xmax = 1.05;

figure();
barh(XAxis, YAxis1);
colororder(newcolors{1})
title('HC Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\d_Bar_1.png'))
figure();
barh(XAxis, YAxis2);
colororder(newcolors{2})
title('Premanifest Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\d_Bar_2.png'))
figure();
barh(XAxis, YAxis3);
colororder(newcolors{3})
title('Mild Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\d_Bar_3.png'))
figure();
barh(XAxis, YAxis4);
colororder(newcolors{4})
title('Severe Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\d_Bar_4.png'))
figure();
barh(XAxis, YAxis5);
colororder(newcolors{5})
title('HGPS Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\d_Bar_5.png'))
%
%$%$ BAR PLOT CAP
myCellTotal = [myCellHCCAPS; myCellRCAPS; myCellMCAPS; myCellSCAPS; myCellPCAPS];
XAxis = categorical(myCellTotal(:,1));
XAxis = reordercats(XAxis,flip(myCellTotal(:,1)));
YAxis = str2double(myCellTotal(:,2:end));
% sumY = sum(YAxis, 2);
% YAxis = bsxfun(@rdivide,YAxis,sumY(:));

Xmax = max(0.5);
for item=1:size(myCellTotal,2)-1
%     YAxis = Answer(:,item);
    
    YAxisPlot = YAxis(:,item)/sum(YAxis(:,item));
%     Xmax = max(YAxisPlot);
    figure();
    barh(XAxis, YAxisPlot);
    colororder(newcolors{6})
    title(uniqueCAPS(item));
    if item == 1
        title('HC');
    end
    if item == size(myCellTotal,2)-1
        title('HGPS');
    end
    xlim([0 Xmax])
    set(gca,'FontSize',14, 'fontweight','bold')
    g = gcf;
    saveas(g, strcat(folderCurrent, '\e_Bar_' + string(item) +'.png'))
end
% %% Clusters Wells
% runInsiders = false;
% if runInsiders
%     folderCurrent = folderCurrent + 'Insiders\';
%     mkdir(folderCurrent);
% end
% 
% temp_CAPPatient = string(DataFileWell.CAPPatient(combinedIndex));
% 
% cMap = [linspace(0,1,128)'*[1,1], ones(128,1)]; % The blue range
% cMap = [cMap; rot90(cMap,2)];                   % add the red range
% cMap = (cMap);
%     
% rowLabels = string(temp_cellType);
% for i=1:size(rowLabels, 1)
%     rowLabels(i) = rowLabels(i) + '_' + string(i) + '_' + temp_CAPPatient(i);
% end
% cg = clustergram(normalized_features, ...
%                              'ColumnLabels', feature_names,...
%                              'RowLabels', rowLabels,...
%                              'RowPdist', 'cityblock',...
%                              'ColumnPdist', 'cityblock',...
%                              'Linkage', 'ward',...
%                              'Cluster', 'COLUMN',...
%                              'Dendrogram', 0,'Colormap', cMap,...
%                            "Standardize", 'column');
% %
% % GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12; gi13; gi14; gi15; gi16];
% % GI = [gi1; gi2; gi5; gi6; gi7; gi8; gi9];
% GINums = [6004 5999 5978 5993 5997 5990 5996];
% % GINums = [6004 5999 5998 6001 5978];
% % GINums = [6004 5975 5976 5989 5954 5980 5984 5992 5978 6001]; % 10
% % GINums = [6004 5975 5976 5989 5954 5980 5984 5992 5993 5983 5974 5979 5978]; %good
% 
% % GINums = [6004 5975 5976 5964 5960 5954 5980 5984 5952 5971 5985 5987 5983 5974 5979 5978];
% 
% %4
% % GI = [gi1; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12; gi13; gi14; gi15];
% % GINums = [4506 4475 4456 4472 4436 4477 4494 4483 4488 4491 4489 4479 4460 4459];
% %3
% % GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12; gi13; gi14; gi15; gi16; gi17];
% % GINums = [4590 4579 4571 4575 4562 4576 4572 4581 4573 4550 4559 4557 4564 4536 4548 4558 4552];
% % GINums = [4590 4579 4571 4575 4562 4576 4572 4581 4573 4550 4559];
% %2
% % GI = [ng1; ng2; ng3; ng4; ng5; ng6; ng7; ng8; ng9; ng10; gi11; gi12; gi13];
% % GINums = [4786 4756 4752 4769 4762 4766 4761 4768 4758 4767 4772 4775 4765];
% %1
% % GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8; gi9; gi10; gi11; gi12];
% % GINums = [5927 5909 5914 5911 5903 5912 5898 5905 5907 5913 5910 5900];
% % GI = [gi1; gi2; gi3; gi4; gi5; gi6; gi7; gi8];
% % GINums = [6484 6455 6470 6475 6477 6474 6476 6473];
% %
% % cluster Group
% de = 1;
% seperate = 0; %33 32 42 41 40
% 
% clear groupsCell
% numGroups = length(cg.RowLabels) - 1;
% numGroups = length(GI);
% clusterDE = {};
% clusterHC = {};
% clusterR = {};
% clusterM = {};
% clusterS = {};
% clusterP = {};
% clusterLabels = {};
% clusterNames = {};
% myCellHC = {};
% myCellR = {};
% myCellM = {};
% myCellS = {};
% myCellP = {};
% myCellHCCAPS = {};
% myCellRCAPS = {};
% myCellMCAPS = {};
% myCellSCAPS = {};
% myCellPCAPS = {};
% 
% % groupsCell = cell([numGroups 1]);
% groupIndexed=0;
% uniqueCAPS = natsort(unique(regexp(string(data.CAPPatient2(combinedIndex)),'\d*','Match', 'once')));
% saved_var = cell(1, numGroups);
% 
% for groupI=numGroups:-1:1
%     group = GINums(groupI);
%     myBreakFlag = false;
%     flag = false;
%     currentGroupInfo = getGroupInfo(cg,group, 1);
%     currentGroupInfo = GI(groupI);
%     groupsCell(group) = currentGroupInfo;
%     if ~runInsiders
%         for iterGroup=1:length(clusterLabels)
%             currNames = groupsCell(clusterLabels{iterGroup}).RowNodeNames;
%             if sum(ismember(currNames, currentGroupInfo.RowNodeNames)) ~= 0
%                 myBreakFlag = true;
%                 break;
%             end
%         end
%         if myBreakFlag
%            continue; 
%         end
%     end
%     countSevere = sum(contains(string(currentGroupInfo.RowNodeNames), "Severe"));
%     countMild = sum(contains(string(currentGroupInfo.RowNodeNames), "Mild"));
%     countRisk = sum(contains(string(currentGroupInfo.RowNodeNames), "Premanifest"));
%     countHC = sum(contains(string(currentGroupInfo.RowNodeNames), "HC"));
%     countProgeria = sum(contains(string(currentGroupInfo.RowNodeNames), "HGPS"));
%     
%     CAPScores = containers.Map;
%     for uC=1:length(uniqueCAPS)
%         CAPScores(uniqueCAPS(uC)) = 0;
%     end
%     inner_values = zeros(1, length(currentGroupInfo.RowNodeNames));
%     for item=1:length(currentGroupInfo.RowNodeNames)
%         currString = string(currentGroupInfo.RowNodeNames(item));
%         pattern = '\d+';
%         matches = regexp(currString, pattern, 'match');
%         first_number = str2double(matches{1}); % Convert the first match to a number
%         inner_values(item) = first_number;
% 
%         split_string = strsplit(currString, '_');
%         resultCAP = regexp(split_string{3},'\d*','Match', 'once');
%         CAPScores(resultCAP) = CAPScores(resultCAP) + 1;
%     end
%     saved_var{groupI} = inner_values;
% 
%     [~,sortedCAPScoresIndexes] = natsort(CAPScores.keys);
%     CAPScoresString = string(CAPScores.values);
%     CAPScoresString = CAPScoresString(sortedCAPScoresIndexes);
%     
%     if countMild < countHC && countRisk < countHC && countSevere < countHC && countProgeria < countHC% HC
%     if max([countMild countRisk countSevere countProgeria]) / countHC < de
%         if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%             flag = true;
%             numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
%             if sum(ismember(numbers, group)) ~= 0
% %                continue; 
%             end
%             clusterLabels(end+1) = {group};
%             clusterDE(end+1) = {max([countMild countRisk countSevere countProgeria]) / countHC};
%             clusterHC(end+1) = {countHC};
%             clusterR(end+1) = {countRisk};
%             clusterM(end+1) = {countMild};
%             clusterS(end+1) = {countSevere};
%             clusterP(end+1) = {countProgeria};
%         end
%     end
%     end
%     if countMild < countSevere && countRisk < countSevere && countHC < countSevere && countProgeria < countSevere% Severe
%     if max([countMild countRisk countHC countProgeria]) / countSevere < de
%         if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%             flag = true;
%             numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
%             if sum(ismember(numbers, group)) ~= 0
% %                continue; 
%             end
%             clusterLabels(end+1) = {group};
%             clusterDE(end+1) = {max([countMild countRisk countHC countProgeria]) / countSevere};
%             clusterHC(end+1) = {countHC};
%             clusterR(end+1) = {countRisk};
%             clusterM(end+1) = {countMild};
%             clusterS(end+1) = {countSevere};
%             clusterP(end+1) = {countProgeria};
%         end
%     end
%     end
% 
%     if countSevere < countMild && countRisk < countMild && countHC < countMild && countProgeria < countMild% Mild
%     if max([countSevere countRisk countHC countProgeria]) / countMild < de
%         if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%             flag = true;
%             numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
%             if sum(ismember(numbers, group)) ~= 0
% %                continue; 
%             end
%             clusterLabels(end+1) = {group};
%             clusterDE(end+1) = {max([countSevere countRisk countHC countProgeria]) / countMild};
%             clusterHC(end+1) = {countHC};
%             clusterR(end+1) = {countRisk};
%             clusterM(end+1) = {countMild};
%             clusterS(end+1) = {countSevere};
%             clusterP(end+1) = {countProgeria};
%         end
%     end
%     end
%     
%     if countSevere < countRisk && countMild < countRisk && countHC < countRisk && countProgeria < countRisk% Risk
%     if max([countSevere countMild countHC countProgeria]) / countRisk < de
%         if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%             flag = true;
%             numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
%             if sum(ismember(numbers, group)) ~= 0
% %                continue; 
%             end
%             clusterLabels(end+1) = {group};
%             clusterDE(end+1) = {max([countSevere countMild countHC countProgeria]) / countRisk};
%             clusterHC(end+1) = {countHC};
%             clusterR(end+1) = {countRisk};
%             clusterM(end+1) = {countMild};
%             clusterS(end+1) = {countSevere};
%             clusterP(end+1) = {countProgeria};
%         end
%     end
%     end
%     
%     if countSevere < countProgeria && countMild < countProgeria && countHC < countProgeria && countRisk < countProgeria% Progeria
%     if max([countSevere countMild countHC countRisk]) / countProgeria < de
%         if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%             flag = true;
%             numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
%             if sum(ismember(numbers, group)) ~= 0
% %                continue; 
%             end
%             clusterLabels(end+1) = {group};
%             clusterDE(end+1) = {max([countSevere countMild countHC countRisk]) / countProgeria};
%             clusterHC(end+1) = {countHC};
%             clusterR(end+1) = {countRisk};
%             clusterM(end+1) = {countMild};
%             clusterS(end+1) = {countSevere};
%             clusterP(end+1) = {countProgeria};
%         end
%     end
%     end
%     
%         if flag
%             disp("******************************************************");
%             checkFlag = false;
%         if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countHC
%             checkFlag = true;
%            disp(strcat(string(group), "_HC_Cluster"));
%            disp(strcat("Num_Of_HC_", string(countHC)));
%            disp(strcat("Num_Of_Severe_", string(countSevere)));
%            disp(strcat("Num_Of_Mild_", string(countMild)));
%            disp(strcat("Num_Of_Premanifest_", string(countRisk)));
%            disp(strcat("Num_Of_HGPS_", string(countProgeria)));
%            if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%                groupIndexed = groupIndexed + 1;
%                 clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-HC")};
%                 myCellHC = [myCellHC; strcat('Cluster-',string(groupIndexed), "-HC"), countHC, countRisk, countMild, countSevere, countProgeria];
%                 myCellHCCAPS = [myCellHCCAPS; strcat('Cluster-',string(groupIndexed), "-HC"), CAPScoresString];
% 
%            end
%         end
%         if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countSevere
%             if checkFlag
%                  disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
%             end
%             checkFlag = true;
%            disp(strcat(string(group), "_Severe_Cluster"));
%            disp(strcat("Num_Of_HC_", string(countHC)));
%            disp(strcat("Num_Of_Severe_", string(countSevere)));
%            disp(strcat("Num_Of_Mild_", string(countMild)));
%            disp(strcat("Num_Of_Premanifest_", string(countRisk)));
%            disp(strcat("Num_Of_HGPS_", string(countProgeria)));
%            if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%                groupIndexed = groupIndexed + 1;
%                 clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-Severe")};
%                 myCellS = [myCellS; strcat('Cluster-',string(groupIndexed), "-Severe"), countHC, countRisk, countMild, countSevere, countProgeria];
%                 myCellSCAPS = [myCellSCAPS; strcat('Cluster-',string(groupIndexed), "-Severe"), CAPScoresString];
%            end
%         end
%         if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countMild
%             if checkFlag
%                  disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
%             end
%             checkFlag = true;
%            disp(strcat(string(group), "_Mild_Cluster"));
%            disp(strcat("Num_Of_HC_", string(countHC)));
%            disp(strcat("Num_Of_Severe_", string(countSevere)));
%            disp(strcat("Num_Of_Mild_", string(countMild)));
%            disp(strcat("Num_Of_Premanifest_", string(countRisk)));
%            disp(strcat("Num_Of_Progeria_", string(countProgeria)));
%            if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%                groupIndexed = groupIndexed + 1;
%                 clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-Mild")};
%                 myCellM = [myCellM; strcat('Cluster-',string(groupIndexed), "-Mild"), countHC, countRisk, countMild, countSevere, countProgeria];
%                 myCellMCAPS = [myCellMCAPS; strcat('Cluster-',string(groupIndexed), "-Mild"), CAPScoresString];
% 
%            end
%         end
%         if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countRisk
%             if checkFlag
%                  disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
%             end
%             checkFlag = true;
%            disp(strcat(string(group), "_Risk_Cluster"));
%            disp(strcat("Num_Of_HC_", string(countHC)));
%            disp(strcat("Num_Of_Severe_", string(countSevere)));
%            disp(strcat("Num_Of_Mild_", string(countMild)));
%            disp(strcat("Num_Of_Premanifest_", string(countRisk)));
%            disp(strcat("Num_Of_HGPS_", string(countProgeria)));
%            if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%                groupIndexed = groupIndexed + 1;
%                 clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-Premanifest")};
%                 myCellR = [myCellR; strcat('Cluster-',string(groupIndexed), "-Premanifest"), countHC, countRisk, countMild, countSevere, countProgeria];
%                 myCellRCAPS = [myCellRCAPS; strcat('Cluster-',string(groupIndexed), "-Premanifest"), CAPScoresString];
%             end
%         end
%         if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countProgeria
%             if checkFlag
%                  disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
%             end
%             checkFlag = true;
%            disp(strcat(string(group), "_HGPS_Cluster"));
%            disp(strcat("Num_Of_HC_", string(countHC)));
%            disp(strcat("Num_Of_Severe_", string(countSevere)));
%            disp(strcat("Num_Of_Mild_", string(countMild)));
%            disp(strcat("Num_Of_Premanifest_", string(countRisk)));
%            disp(strcat("Num_Of_HGPS_", string(countProgeria)));
%            if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
%                groupIndexed = groupIndexed + 1;
%                 clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-HGPS")};
%                 myCellP = [myCellP; strcat('Cluster-',string(groupIndexed), "-HGPS"), countHC, countRisk, countMild, countSevere, countProgeria];
%                 myCellPCAPS = [myCellPCAPS; strcat('Cluster-',string(groupIndexed), "-HGPS"), CAPScoresString];
%             end
%         end
%         
%         disp(strcat("Total_Cells_Is_", string(countMild + countRisk + countSevere + countHC + countProgeria)));
% 
%         end
% end
%     
% groupsCellMotility = groupsCell;
% clusterNamesMotility = clusterNames;
% %
% CM = jet(length(clusterLabels));
% CM(:,4) = 0.5;
% CMCell = {};
% for color=1:size(CM,1)
%     CMCell(end+1) = {CM(color,:)};
% end
% CMCell = CMCell(randperm(numel(CMCell)));
% CMCell = {
%     [0, 0.2, 0.4],...
%     [0.3, 0.8, 1],...
%     [1, 0.4, 0.4],...
%     [0.8500, 0.3250, 0.0980],...
%     [0.5, 0, 0],...
%     [0, 0.4470, 0.7410] ,...
%     [128, 0, 128]/255,...
% };
% %
% if ~isempty(clusterLabels)
% rm = struct('GroupNumber',clusterLabels,'Annotation',clusterNames,...
%     'Color',CMCell, 'FontSize', 5);
% set(cg,'RowGroupMarker',rm)
% cgf = plot(cg);
% set(cgf,'FontSize',14, 'fontweight','bold')
% fig = get(cgf, 'Parent');
% fig.WindowState = 'maximized';
% saveas(cgf,strcat(folderCurrent, '\c_ClustergramWithClusterNames.png'))
% end
%
close all force
%$%$ Draw
if ~isempty(clusterLabels)
    
    for iBig=1:1:length(clusterLabels)
    starts = iBig;
    ends = iBig+0;
    if ends > length(clusterLabels)
        ends = length(clusterLabels); 
    end
    fig = figure('Position', get(0, 'Screensize'));
    t = tiledlayout(3,1, 'TileSpacing','tight');
    for myGroup=starts:ends
        currentGroupName = clusterNames{myGroup};
        if strfind(currentGroupName, 'Cluster 6') > 0
            myColorO = newcolors{2};
            myTypeOuter = 'Risk';
        elseif strfind(currentGroupName, 'Mild') > 0
            myColorO = newcolors{3};
            myTypeOuter = 'Mild';
        elseif (contains(currentGroupName, 'Cluster 2') || contains(currentGroupName, 'Cluster 3') || contains(currentGroupName, 'Cluster 4'))
            myColorO = newcolors{4};
            myTypeOuter = 'Severe';
        elseif strfind(currentGroupName, 'Cluster 1') > 0
            myColorO = newcolors{1};
            myTypeOuter = 'HC';
        elseif (contains(currentGroupName, 'Cluster 5') || contains(currentGroupName, 'Cluster 7'))
            myColorO = newcolors{5};
            myTypeOuter = 'HGPS';
        end
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        aPol = [];
        aColors = [];
        aType = [];
        RiskCount = 0;
        MildCount = 0;
        SevereCount = 0;
        HCCount = 0;
        ProgeriaCount = 0;
        for lookedGroup=1:length(currentLookedGroup)
            currentWell = string(currentLookedGroup(lookedGroup));
            if strfind(currentWell, 'Premanifest') > 0
                RiskCount = RiskCount + 1;
                myColor = newcolors{2};
                myType = 'Risk';
            elseif strfind(currentWell, 'Mild') > 0
                MildCount = MildCount + 1;
                myColor = newcolors{3};
                myType = 'Mild';
            elseif strfind(currentWell, 'Severe') > 0
                SevereCount = SevereCount + 1;
                myColor = newcolors{4};
                myType = 'Severe';
            elseif strfind(currentWell, 'HC') > 0
                HCCount = HCCount + 1;
                myColor = newcolors{1};
                myType = 'HC';
            elseif strfind(currentWell, 'HGPS') > 0
                ProgeriaCount = ProgeriaCount + 1;
                myColor = newcolors{5};
                myType = 'HGPS';
            end
            number = str2double(regexp(currentWell,'\d*','Match', 'once'));
            aPol = [aPol; xys(number)];
            vColors = repmat({myColor},size(xys(number),1),1);
            aColors = [aColors; vColors];
            vType = repmat({myType},size(xys(number),1),1);
            aType = [aType; vType];
        end
        if size(aPol,1) < numOfRandoms
            numOfRandomsIn = size(aPol,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
        out = [];
        indexBuild = 0;
        typeOuter = 0;
        typeElse = 0;
        if strcmp(myTypeOuter, 'Risk')
            usedOuter = RiskCount;
            needed = RiskCount / max([MildCount SevereCount HCCount ProgeriaCount]);
        elseif strcmp(myTypeOuter, 'Mild')
            usedOuter = MildCount;
            needed = MildCount / max([RiskCount SevereCount HCCount ProgeriaCount]);
        elseif strcmp(myTypeOuter, 'Severe')
            usedOuter = SevereCount;
            needed = SevereCount / max([MildCount RiskCount HCCount ProgeriaCount]);
        elseif strcmp(myTypeOuter, 'HC')
            usedOuter = HCCount;
            needed = HCCount / max([MildCount SevereCount RiskCount ProgeriaCount]);
        elseif strcmp(myTypeOuter, 'HGPS')
            usedOuter = ProgeriaCount;
            needed = ProgeriaCount / max([MildCount SevereCount RiskCount HCCount]);
        end
        calc = 1 - 1/needed;
        breakLoop = false;
        while true
            if breakLoop
                break;
            end
            if size(out, 2) == numOfRandomsIn
                breakLoop = true;
                break;
            end
            r = randi([1 size(aPol,1)],1,1);
            if sum(ismember(out, r)) ~= 0
                continue;
            end
            if strcmp(aType(r),myTypeOuter)
               
                if typeOuter > ceil(numOfRandomsIn*calc) + 2
                    continue;
                end
                typeOuter = typeOuter + 1;
                out = [out r];
            else
                if typeOuter > ceil(numOfRandomsIn*calc) + 2 || typeOuter == usedOuter
                   typeElse = typeElse + 1;
                    out = [out r]; 
                end
            end
            
            
        end
        out = out(randperm(length(out)));
%         out = randperm(size(aPol,1),numOfRandomsIn);
        outAll = randperm(size(aPol,1),size(aPol,1));
        counter_x = 0;
        counter_y = 0;
        space = 500;
        nexttile
        hold on
        for j=1:numOfRandomsIn
            currNum = out(j);
            % getting sepcific cell trajectory
            xy=aPol{currNum}; 
            col = aColors{currNum};

            % zero cell position at first frame.
            xy=xy-ones(size(xy(:,1)))*xy(1,:);
            plot(xy(:,1)+counter_x*space,xy(:,2)+counter_y*space,'Color',cell2mat(CMCell(myGroup)), 'LineWidth', 2);
            counter_x = counter_x + 1;
            if mod(counter_x,10)==0
                counter_x = 0;
                counter_y = counter_y + 1;
            end
        end
        set(gca,'xtick',[],'ytick',[]);
        title({"Trajectories", string(clusterNames{myGroup})});
        axis square;
        hold off
        box on
        ylim([-space space*5])
        xlim([-space space*10])
    % Show Sunplot of Cells
        nexttile
        cidk = jet(numOfRandomsIn);
        for k=1:1:numOfRandomsIn
             currNum = out(k);
             xy=aPol{currNum};         
             xy=xy-ones(size(xy(:,1)))*xy(1,:); % zero cell position at first frame.
              xy1 = smooth(xy(:,1));
             xy2 = smooth(xy(:,2));
             plot(xy1,xy2,'-','color',[cidk(k,:) 1],'linewidth', 2); hold on; 
        end
%         title(strcat("Sunplot Trajectories for ", string(clusterNames{myGroup})));
        box on
        axis square
        hold off;
        set(gca,'xtick',[],'ytick',[]);
        ylim([-350 350])
        xlim([-350 350])
        text(90, 300, strcat("N = ", string(length(aPol))),'FontSize', 16, 'FontWeight', 'bold');

        
        % Show merged of Cells
        nexttile
        cidk = jet(numOfRandomsIn);
        xy = 0;
        for k=1:1:numOfRandomsIn
             currNum = out(k);
             xyC=aPol{currNum};         
             xy=xy + xyC-ones(size(xyC(:,1)))*xyC(1,:); % zero cell position at first frame.
        end
        xy = xy / numOfRandomsIn;
        xy1 = smooth(xy(:,1));
        xy2 = smooth(xy(:,2));
        plot(xy1,xy2,'-','color',cell2mat(CMCell(myGroup)), 'linewidth',2);
%         title(strcat("Merged Trajectories for ", string(clusterNames{myGroup})));
        box on
        axis square
        hold off;
        set(gca,'xtick',[],'ytick',[]);
        ylim([-30 30])
        xlim([-30 30])

    end
    exportgraphics(t,strcat(folderCurrent, 'e_', string(myGroup),'_TRJ.jpg'),'Resolution',500)
    end
end
close all force

%$%$ Features
if ~isempty(clusterLabels)

    vNumbersAll = [];
    normedAll = [];
    vMeansAll = [];
    importanceAll = [];
    importanceNormed = [];
for myGroup=1:1:length(clusterLabels)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        vNumbersAll = [vNumbersAll; currNums];
        vMeansAll = [vMeansAll; mean(z(currNums,:))];
end
zClusters = z(vNumbersAll, :);
normedzClusters = rescale(zClusters, 'InputMin',min(zClusters),'InputMax',max(zClusters));

    for iBig=1:1:length(clusterLabels)
    starts = iBig;
    ends = iBig+0;
    if ends > length(clusterLabels)
        ends = length(clusterLabels); 
    end
for myGroup=starts:ends

        currentGroupName = clusterNames{myGroup};
        if strfind(currentGroupName, 'Premanifest') > 0
            myColorO = newcolors{2};
            myTypeOuter = 'Risk';
        elseif strfind(currentGroupName, 'Mild') > 0
            myColorO = newcolors{3};
            myTypeOuter = 'Mild';
        elseif strfind(currentGroupName, 'Severe') > 0
            myColorO = newcolors{4};
            myTypeOuter = 'Severe';
        elseif strfind(currentGroupName, 'HC') > 0
            myColorO = newcolors{1};
            myTypeOuter = 'HC';
        elseif strfind(currentGroupName, 'HGPS') > 0
            myColorO = newcolors{5};
            myTypeOuter = 'HGPS';
        end
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        vNumbers = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(vNumbers) = 1;
%         aPol = [];
%         aColors = [];
%         aType = [];
        aFeature = [];
        RiskCount = 0;
        MildCount = 0;
        SevereCount = 0;
        HCCount = 0;
    
    % which features are the strongest
    feature_importance = [];
    temp_model = fitcensemble(z,aNumbers, 'Method','bag','NumLearningCycles', 1000);
    importance = oobPermutedPredictorImportance(temp_model);
    importanceAll = [importanceAll; importance];
    importanceNormed = [importanceNormed; rescale(importance, 'InputMin',min(importance),'InputMax',max(importance))];
    importance = rescale(importance, 'InputMin',min(importance),'InputMax',max(importance));
    [idxProb,idx] = sort(importance,'descend');
    feature_names_sorted_by_importance = feature_names(idx);
%     feature_importance = [feature_importance; feature_names_sorted_by_importance];
%     feature_importance_top = feature_importance(1:4);
    feature_importance_top = feature_names_sorted_by_importance(idxProb > 0.5);
    fig = figure('Position', get(0, 'Screensize'));
    t = tiledlayout(ceil(size(feature_importance_top,2)/2),2, 'TileSpacing','Compact');   
    %

        
%     for featureToCheck=1:size(feature_names,2)
%         if ~ismember(string(feature_importance_top), string(feature_names(featureToCheck)))
%             continue;
%         end
    featNames = [];
    aFeatureDataAll = [];
    for featureToCheckTop=1:size(feature_importance_top,2)
        featureToCheck = ismember(feature_names, feature_importance_top(featureToCheckTop));
        featureToCheck = find(featureToCheck);
        featNames = [featNames; string(feature_importance_top(featureToCheckTop))];
        aFeature = [];
        for lookedGroup=1:length(currentLookedGroup)
            currentWell = string(currentLookedGroup(lookedGroup));
            if strfind(currentWell, 'Risk') > 0
                RiskCount = RiskCount + 1;
                myColor = newcolors{2};
                myType = 'Risk';
            elseif strfind(currentWell, 'Mild') > 0
                MildCount = MildCount + 1;
                myColor = newcolors{3};
                myType = 'Mild';
            elseif strfind(currentWell, 'Severe') > 0
                SevereCount = SevereCount + 1;
                myColor = newcolors{4};
                myType = 'Severe';
            elseif strfind(currentWell, 'HC') > 0
                HCCount = HCCount + 1;
                myColor = newcolors{1};
                myType = 'HC';
            elseif strfind(currentWell, 'HGPS') > 0
                ProgeriaCount = ProgeriaCount + 1;
                myColor = newcolors{5};
                myType = 'HGPS';
            end
            number = str2double(regexp(currentWell,'\d*','Match', 'once'));
            aFeature = [aFeature; z(number, featureToCheck)];
%             aPol = [aPol; xys(number)];
%             vColors = repmat({myColor},size(xys(number),1),1);
%             aColors = [aColors; vColors];
%             vType = repmat({myType},size(xys(number),1),1);
%             aType = [aType; vType];
        end
        nexttile;
        [f1,xi1] = ksdensity(aFeature, 'Kernel','normal', 'Support',[floor(min(aFeature)) ceil(max(aFeature))], 'NumPoints',size(aFeature,1));
        area(xi1,f1, 'FaceColor', cell2mat(CMCell(myGroup)), 'FaceAlpha', 0.5)
        titleSMR = strcat("Feature Importance", {' '}, num2str(featureToCheckTop),{' '}, "for Feature", {' '}, feature_names(featureToCheck), {' '}, "for", {' '}, string(clusterNames{myGroup}));
        title(titleSMR);
        aFeatureData = [min(aFeature) max(aFeature) mean(aFeature) median(aFeature) std(aFeature)];
        aFeatureDataAll = [aFeatureDataAll; aFeatureData];
        
        

    end
    exportgraphics(t,strcat(folderCurrent, 'f_', string(myGroup),'_FeatureImportance.jpg'),'Resolution',500)
    f = figure();
    h = uitable(f, 'data', aFeatureDataAll, 'ColumnName', {'Min', 'Max', 'Mean', 'Median', 'STD'}, 'RowName', featNames);
%     h.Position = [540 360 700 300];
    h.Position = [40 60 500 300];
%     txt_title = uicontrol('Style', 'text', 'Position', [540 460 700 300], 'String', string(clusterNames{myGroup}));
    saveas(h,strcat(folderCurrent, 'i_', string(myGroup),'_Table.jpg'))

    % Radar Total
    normedZ = rescale(z, 'InputMin',min(z),'InputMax',max(z));
    r = figure;
    normed = normedZ(vNumbers, :);
%     normed = rescale(normed, 'InputMin',min(z),'InputMax',max(z));
    normed = mean(normed);
    zerosS = 0*ones(1,size(normed,2));
    onesS = ones(1,size(normed,2));
    L = [zerosS;onesS];
    spider_plot(normed,'AxesLimits',L, ...
        'AxesPrecision', onesS,...
        'FillOption', {'on'},...
        'Color', cell2mat(CMCell(myGroup)), ...
        'FillTransparency', 0.05,...
        'AxesLabels', feature_names,...
        'AxesLabelsColors', false(1,9), ...
        'AxesInterval', 5);
%         'AxesLabelsColors', ismember(string(feature_names), string(feature_importance_top))
    title(strcat('Spider Plot (by All Data) for',{' '}, string(clusterNames{myGroup})));
    exportgraphics(r,strcat(folderCurrent, 'g_', string(myGroup),'_RadarTotal.jpg'),'Resolution',500)
    
%     myText = sprintf("I have a cluster with %f Dp, %f Dtot, %f Pp, %f Psi, %f Pnp, %f Dnp, %f MSD, %f Sp, %f Snp. What could it mean?", normed(1), normed(2), normed(3), normed(4), normed(5), normed(6), normed(7), normed(8), normed(9));
%     disp(myText);
    %
    % Radar Local
    r = figure;
    normed = z(vNumbers, :);
    normed = mean(normed);
    normed = rescale(normed, 'InputMin',min(vMeansAll),'InputMax',max(vMeansAll));
%     normedAll = [normedAll; normed];
    zerosS = 0*ones(1,size(normed,2));
    onesS = ones(1,size(normed,2));
    L = [zerosS;onesS];
    spider_plot(normed,'AxesLimits',L, ...
        'AxesPrecision', onesS,...
        'FillOption', {'on'},...
        'Color', cell2mat(CMCell(myGroup)), ...
        'FillTransparency', 0.05,...
        'AxesLabels', feature_names,...
        'AxesLabelsColors', false(1,9), ...
        'AxesInterval', 5);
    %         'AxesLabelsColors', ismember(string(feature_names), string(feature_importance_top))

    title(strcat('Spider Plot (by Clusters Data) for',{' '}, string(clusterNames{myGroup})));
    exportgraphics(r,strcat(folderCurrent, 'h_', string(myGroup),'_RadarClusters.jpg'),'Resolution',500)
    
    myText = sprintf("I have a cluster with %f Dp, %f Dtot, %f Pp, %f Psi, %f Pnp, %f Dnp, %f MSD, %f Sp, %f Snp. What could it mean?", normed(1), normed(2), normed(3), normed(4), normed(5), normed(6), normed(7), normed(8), normed(9));
    disp(myText);
end
    end
end
heatmapStr = ["Total Cells", "HC Cells", "Premanifest Cells", "Mild Cells", "Severe Cells", "HGPS Cells", "Seperation Factor"];
% figure();
allClusters = [cell2mat(clusterS)+cell2mat(clusterR)+cell2mat(clusterM)+cell2mat(clusterHC)+cell2mat(clusterP); cell2mat(clusterHC); cell2mat(clusterR);cell2mat(clusterM);cell2mat(clusterS);cell2mat(clusterP);1-cell2mat(clusterDE)];
% heatmap(clusterNames, heatmapStr, allClusters, 'ColorScaling', 'scaledrows', 'Colormap', flipud(hot), 'ColorbarVisible', 'off')
fig = figure('Position', get(0, 'Screensize'));
t = tiledlayout(2,1, 'TileSpacing','Compact');
% figure();
nexttile;
heatmap(clusterNames, feature_names, importanceAll');
title('Heatmap of feature importance by all data')
nexttile;
heatmap(clusterNames, heatmapStr, allClusters, 'ColorScaling', 'scaledrows', 'Colormap', flipud(hot), 'ColorbarVisible', 'off')
title('Heatmap of data distribuation')
% g = gcf;
% t.WindowState = 'maximized';
saveas(fig,strcat(folderCurrent, 'j_HeatmapAll.jpg'))

fig = figure('Position', get(0, 'Screensize'));
t = tiledlayout(2,1, 'TileSpacing','Compact');
% figure();
nexttile
heatmap(clusterNames, feature_names, importanceNormed');
title('Heatmap of feature importance by Clusters')
% g = gcf;
nexttile
heatmap(clusterNames, heatmapStr, allClusters, 'ColorScaling', 'scaledrows', 'Colormap', flipud(hot), 'ColorbarVisible', 'off')
title('Heatmap of data distribuation')
% t.WindowState = 'maximized';
saveas(fig,strcat(folderCurrent, 'j_HeatmapNormed.jpg'))
close all
% Plot Features {'Dp','Dtot','Pp','Psi','Pnp','Dnp','MSD10','Sp','Snp'}

%% Prediction Intervals
temp_cellType = DataFileWell.CellType(combinedIndex);
patients = DataFileWell.CAPPatient(combinedIndex);
normalized_features = log(feature_vectors);
z = zscore(normalized_features);

featuresToUse = z;

featuresHC = featuresToUse(IDX == 1, :);
patientsHC = patients(IDX == 1);
patientsHCUnique = unique(patientsHC);
PfeaturesHC = zeros(size(patientsHCUnique, 1),9);
for i = 1:size(patientsHCUnique, 1)
   meanPatient = mean(featuresHC(patientsHC == patientsHCUnique(i), :));
   PfeaturesHC(i,:)= meanPatient;
end
featuresSevere = featuresToUse(IDX == 2, :);
patientsSevere = patients(IDX == 2);
patientsSevereUnique = unique(patientsSevere);
PfeaturesSevere = zeros(size(patientsSevereUnique, 1),9);
for i = 1:size(patientsSevereUnique, 1)
   meanPatient = mean(featuresSevere(patientsSevere == patientsSevereUnique(i), :));
   PfeaturesSevere(i,:)= meanPatient;
end
featuresMild = featuresToUse(IDX == 3, :);
patientsMild = patients(IDX == 3);
patientsMildUnique = unique(patientsMild);
PfeaturesMild = zeros(size(patientsMildUnique, 1),9);
for i = 1:size(patientsMildUnique, 1)
   meanPatient = mean(featuresMild(patientsMild == patientsMildUnique(i), :));
   PfeaturesMild(i,:)= meanPatient;
end
featuresPre = featuresToUse(IDX == 4, :);
patientsPre = patients(IDX == 4);
patientsPreUnique = unique(patientsPre);
PfeaturesPre = zeros(size(patientsPreUnique, 1),9);
for i = 1:size(patientsPreUnique, 1)
   meanPatient = mean(featuresPre(patientsPre == patientsPreUnique(i), :));
   PfeaturesPre(i,:)= meanPatient;
end
featuresHGPS = featuresToUse(IDX == 5, :);
patientsHGPS = patients(IDX == 5);
patientsHGPSUnique = unique(patientsHGPS);
PfeaturesHGPS = zeros(size(patientsHGPSUnique, 1),9);
for i = 1:size(patientsHGPSUnique, 1)
   meanPatient = mean(featuresHGPS(patientsHGPS == patientsHGPSUnique(i), :));
   PfeaturesHGPS(i,:)= meanPatient;
end
featuresHGPS2 = featuresToUse(IDX == 6, :);
patientsHGPS2 = patients(IDX == 6);
patientsHGPSUnique2 = unique(patientsHGPS2);
PfeaturesHGPS2 = zeros(size(patientsHGPSUnique2, 1),9);
for i = 1:size(patientsHGPSUnique2, 1)
   meanPatient = mean(featuresHGPS2(patientsHGPS2 == patientsHGPSUnique2(i), :));
   PfeaturesHGPS2(i,:)= meanPatient;
end
featuresHGPS3 = featuresToUse(IDX == 7, :);
patientsHGPS3 = patients(IDX == 7);
patientsHGPSUnique3 = unique(patientsHGPS3);
PfeaturesHGPS3 = zeros(size(patientsHGPSUnique3, 1),9);
for i = 1:size(patientsHGPSUnique3, 1)
   meanPatient = mean(featuresHGPS3(patientsHGPS3 == patientsHGPSUnique3(i), :));
   PfeaturesHGPS3(i,:)= meanPatient;
end

newcolorsRGB = {
    [0, 0.2, 0.4],...
    [0.3, 0.8, 1],...
    [1, 0.4, 0.4],...
    [0.8500, 0.3250, 0.0980],...
    [0.5, 0, 0],...
    [0, 0.4470, 0.7410] ,...
    [128, 0, 128]/255,...
};
confidenceZ = 1.5;

for feature=1:size(featuresToUse, 2)
    figure;
    hold on;
    
    % create a cell array to store features of each group
%     featuresCell = {featuresHC, featuresPre, featuresMild, featuresSevere, featuresHGPS, featuresHGPS2, featuresHGPS3};
    featuresCell = {PfeaturesHC, PfeaturesPre, PfeaturesMild, PfeaturesSevere, PfeaturesHGPS, PfeaturesHGPS2, PfeaturesHGPS3};

    % create a cell array to store group names
    groupName = {'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7'};
    h = cell(1,7);

    for i = 1:7 % for each group
        disp("**********************");
        disp(groupName{i});
        
        featureUse = featuresCell{i};
        
        disp("Feature is:");
        disp(string(feature_names(feature)));
        
        currentFeatureMean = mean(featureUse(:, feature));
        currentFeatureSTD = std(featureUse(:, feature));
        
        disp("Lower Bound is:");
        lowerBound = currentFeatureMean - confidenceZ*currentFeatureSTD;
        disp(string(lowerBound));
        
        disp("Upper Bound is:");
        upperBound = currentFeatureMean + confidenceZ*currentFeatureSTD;
        disp(string(upperBound));
        
        % plot the histogram (you may want to adjust the number of bins)
        [N,edges] = histcounts(featureUse(:, feature), 'Normalization', 'probability');
        edges = edges(1:end-1) + diff(edges)/2;
        h{i} = plot(edges, N, 'Color', newcolorsRGB{i});
        
        % plot the lower and upper bounds
        line([lowerBound lowerBound], get(gca, 'YLim'), 'Color', newcolorsRGB{i}, 'LineStyle', '--');
        line([upperBound upperBound], get(gca, 'YLim'), 'Color', newcolorsRGB{i}, 'LineStyle', '--');
    end
    
    title(string(feature_names(feature)))
    xlabel('Feature Value')
    ylabel('Probability')
    legend([h{1} h{2} h{3} h{4} h{5} h{6} h{7}], groupName)
    
    hold off;
end

%% Z-Factor
% temp_cellType = DataFileWell.CellType(combinedIndex);
% normalized_features = log(feature_vectors);
% z = zscore(normalized_features);

% featuresToUse = z;


% featuresHC = featuresToUse(temp_cellType == "'HC'", :);
% featuresSevere = featuresToUse(temp_cellType == "'Severe'", :);
% featuresMild = featuresToUse(temp_cellType == "'Mild'", :);
% featuresPre = featuresToUse(temp_cellType == "'Premanifest'", :);
% featuresHGPS = featuresToUse(temp_cellType == "'HGPS'", :);

for feature=1:size(featuresToUse, 2)
    
    % create a cell array to store features of each group
%     featuresCell = {featuresHC, featuresPre, featuresMild, featuresSevere, featuresHGPS, featuresHGPS2, featuresHGPS3};
    
    % create a cell array to store group names
%     groupName = {'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7'};
    h = cell(1,7);

    for i = 1:7 % for each group
        disp("**********************");
        disp("Group Outside is:");
        disp(groupName{i});
        
        featureUse = featuresCell{i};
        
        disp("Feature is:");
        disp(string(feature_names(feature)));
        
        currentFeatureMean = mean(featureUse(:, feature));
        currentFeatureSTD = std(featureUse(:, feature));
        for j = 1:5 % for each group
            if j == i
               continue; 
            end
            featureUseInside = featuresCell{j};

            disp("Group Inside is:");
            disp(groupName{j});

            insideFeatureMean = mean(featureUseInside(:, feature));
            insideFeatureSTD = std(featureUseInside(:, feature));
            
%             SeperationBand = (insideFeatureMean - 3*insideFeatureSTD) - (currentFeatureMean + 3*currentFeatureSTD);
            SeperationBand = 1.5*(insideFeatureSTD + currentFeatureSTD);
            DynamicRange = abs(insideFeatureMean - currentFeatureMean);
            zFactor = 1 - SeperationBand/DynamicRange;
            
            if zFactor > 0.5
                disp("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
            end
            ("Z-Factor is:");
            disp(zFactor);
        end
    end
end
%% Gini-Simpson per CAP
close all
tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

values = zeros(1, numel(uniqueTempCAP));
dict = containers.Map(uniqueTempCAP, values);

for myGroup=1:1:length(clusterLabels)
    currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
    currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
    aNumbers = zeros(size(z, 1), 1);
    aNumbers(currNums) = 1;
    types = tempCAP(logical(aNumbers));
    
    numericTempCAP = str2double(string(types));
    roundedNumericTempCAP = round(numericTempCAP);
    types = string(roundedNumericTempCAP);
    
    num_all = length(types);
    uniqueTypes = unique(types);
    for type=1:length(uniqueTypes)
        group_size = length(tempCAPCached(tempCAPCached == uniqueTypes(type)));
        num = sum(types == uniqueTypes(type));
        p = num/num_all;
        name = uniqueTypes(type);
        dict(name) = dict(name) + p^2;
    end
end

% Calculate Gini-Simpson index
for key = keys(dict)
    dict(key{1}) = 1 - dict(key{1});
end

% Define the order and colors
ordered_keys = uniqueTempCAP;

newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end

% Create the bar plot with different colors
figure;
hold on;
avg_values = [];
for i = 1:numel(values)
    if i == 1
      col = colors{1};
      avg_values(end+1) = values(i);
    elseif i < 8 % 6
        if i == 2
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{2};
    elseif i < 16 % 8
        if i == 8
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{3};
    elseif i < 21 % 5
        if i == 16
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{4};
    elseif i == 21
        avg_values(end+1) = values(i);
        col = colors{5};
    end
    bar_handle = bar(i, values(i), 'FaceColor', col);
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

avg_values_cached = avg_values;
figure;
hold on;
for i = 1:numel(avg_values)
    if i == 2
        avg_values(i) = avg_values(i) / 6;
    elseif i == 3
        avg_values(i) = avg_values(i) / 8;
    elseif i ==4
        avg_values(i) = avg_values(i) / 5;
    end
    bar_handle = bar(i, avg_values(i), 'FaceColor', colors{i});
end
hold off;

ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};

% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

avg_values = avg_values_cached;
figure;
hold on;
for i = 1:numel(avg_values)
    if i==1 || i ==5
        continue;
    end
    if i == 2
        avg_values(i) = (avg_values(i) + avg_values(1)) / 7;
    elseif i == 3
        avg_values(i) = avg_values(i) / 8;
    elseif i ==4
        avg_values(i) = (avg_values(i) + avg_values(5)) / 6;
    end
    bar_handle = bar(i-1, avg_values(i), 'FaceColor', colors{i});
end
hold off;

ordered_keys = {'HCPremanifest', 'Mild', 'SevereHGPS'};

% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
%% Simpson per CAP
close all
tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

values = zeros(1, numel(uniqueTempCAP));
dict = containers.Map(uniqueTempCAP, values);

for myGroup = 1:1:length(clusterLabels)
    currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
    currNums = str2double(regexp(currentLookedGroup, '\d*', 'Match', 'once'));
    aNumbers = zeros(size(z, 1), 1);
    aNumbers(currNums) = 1;
    types = tempCAP(logical(aNumbers));

    numericTempCAP = str2double(string(types));
    roundedNumericTempCAP = round(numericTempCAP);
    types = string(roundedNumericTempCAP);

    num_all = length(types);
    uniqueTypes = unique(types);
    for type = 1:length(uniqueTypes)
        num_all_new = sum(tempCAPCached == uniqueTypes(type));
        group_size = length(tempCAPCached(tempCAPCached == uniqueTypes(type)));
        num = sum(types == uniqueTypes(type));
        p = num / num_all;
        name = uniqueTypes(type);
        dict(name) = dict(name) + p^2;
    end
end

% Calculate Simpson index
for key = keys(dict)
    dict(key{1}) = 1 - dict(key{1});
end
% Define the order and colors
ordered_keys = uniqueTempCAP;

newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
% P Value
group_entropies = cell(1, length(values));
for i = 1:numel(values)
    group_entropies{i} = values(i);
end
all_entropies = vertcat(group_entropies{:});
% group_indices = arrayfun(@(x) x * ones(size(group_entropies{x})), 1:length(values), 'UniformOutput', false);
% grouping_var = vertcat(group_indices{:});
% [p_value, ~, ~] = kruskalwallis(all_entropies, grouping_var, 'off');
% num_groups = length(all_entropies);
% significant_pairs = [];
% % Loop through all possible pairs of groups
% for i = 1:num_groups
%     for j = i+1:num_groups
%         group1 = all_entropies(i);
%         group2 = all_entropies(j);
%         
%         % Perform the Mann-Whitney U test
%         [p_value, ~] = ranksum(group1, group2);
%         
%         % Check if the p-value is less than 0.05
%         if p_value < 0.05
%             significant_pairs = [significant_pairs; i, j];
%         end
%     end
% end
% % Display the significant pairs
% disp('Significant pairs (p < 0.05):');
% disp(significant_pairs);
% Create the bar plot with different colors
figure;
hold on;
avg_values = [];
for i = 1:numel(values)
    if i == 1
      col = colors{1};
      avg_values(end+1) = values(i);
    elseif i < 8 % 6
        if i == 2
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{2};
    elseif i < 16 % 8
        if i == 8
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{3};
    elseif i < 21 % 5
        if i == 16
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{4};
    elseif i == 21
        avg_values(end+1) = values(i);
        col = colors{5};
    end
    bar_handle = bar(i, values(i), 'FaceColor', col);
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
ordered_keys(1) = 'HC';
ordered_keys(end) = 'HGPS';
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

avg_values_cached = avg_values;
figure;
hold on;
for i = 1:numel(avg_values)
    if i == 2
        avg_values(i) = avg_values(i) / 6;
    elseif i == 3
        avg_values(i) = avg_values(i) / 8;
    elseif i ==4
        avg_values(i) = avg_values(i) / 5;
    end
    bar_handle = bar(i, avg_values(i), 'FaceColor', colors{i});
end
hold off;

ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};

% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

avg_values = avg_values_cached;
figure;
hold on;
for i = 1:numel(avg_values)
    if i==1 || i ==5
        continue;
    end
    if i == 2
        avg_values(i) = (avg_values(i) + avg_values(1)) / 7;
    elseif i == 3
        avg_values(i) = avg_values(i) / 8;
    elseif i ==4
        avg_values(i) = (avg_values(i) + avg_values(5)) / 6;
    end
    bar_handle = bar(i-1, avg_values(i), 'FaceColor', colors{i});
end
hold off;

ordered_keys = {'HCPremanifest', 'Mild', 'SevereHGPS'};

% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

new_group_definitions = {
    [1, 2, 3, 4, 5, 6, 7],  ...
    [8, 9, 10, 11, 12, 13, 14, 15],   ...
    [16, 17, 18, 19, 20, 21]  
};
num_new_groups = length(new_group_definitions);
combined_entropies = cell(num_new_groups, 1);
for i = 1:num_new_groups
    original_group_indices = new_group_definitions{i};
    combined_entropies{i} = all_entropies(original_group_indices);
end
significant_pairs = [];
p_values = [];
for i = 1:num_new_groups
    for j = i+1:num_new_groups
        group1 = combined_entropies{i};
        group2 = combined_entropies{j};
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1, group2);
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            p_values = [p_values;p_value];
        end
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% Shannon per CAP
close all
tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

values = zeros(1, numel(uniqueTempCAP));
dict = containers.Map(uniqueTempCAP, values);
dictCached = containers.Map(uniqueTempCAP, values);


for myGroup=1:1:length(clusterLabels)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        
        num_all = length(types);
        uniqueTypes = unique(types);
        for type=1:length(uniqueTypes)
            num_all_new = sum(tempCAPCached == uniqueTypes(type));
            group_size = length(tempCAPCached(tempCAPCached == uniqueTypes(type)));
            num = sum(types == uniqueTypes(type));
            p = num/num_all;
            name = uniqueTypes(type);
%             name = str2double(name);
%             name = round(name);
%             name = string(name);
%             dict(name) = dict(name) - (p * log(p))/log(group_size);
if type == length(uniqueTypes)
    a = 1;
end
            dict(name) = dict(name) - (p * log(p));
            dictCached(name) = dictCached(name) - (p * log(p));
        end
end
% Define the order and colors
ordered_keys = uniqueTempCAP;

newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
    valuesCached(i) = dictCached(key);
end

% Create the bar plot with different colors
figure;
hold on;
avg_values = [];
for i = 1:numel(values)
    if i == 1
      col = colors{1};
      avg_values(end+1) = values(i);
    elseif i < 8 % 6
        if i == 2
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{2};
    elseif i < 16 % 8
        if i == 8
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{3};
    elseif i < 21 % 5
        if i == 16
            avg_values(end+1) = values(i);
        else
          avg_values(end) = avg_values(end) + values(i);
        end
        col = colors{4};
    elseif i == 21
        avg_values(end+1) = values(i);
        col = colors{5};
    end
    bar_handle = bar(i, values(i), 'FaceColor', col);
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
ordered_keys(1) = 'HC';
ordered_keys(end) = 'HGPS';

xticklabels(ordered_keys);

% P Value
group_entropies = cell(1, length(values));
for i = 1:numel(values)
    group_entropies{i} = values(i);
end
all_entropies = vertcat(group_entropies{:});
% group_indices = arrayfun(@(x) x * ones(size(group_entropies{x})), 1:length(values), 'UniformOutput', false);
% grouping_var = vertcat(group_indices{:});
% [p_value, ~, ~] = kruskalwallis(all_entropies, grouping_var, 'off');
% num_groups = length(all_entropies);
% significant_pairs = [];
% % Loop through all possible pairs of groups
% for i = 1:num_groups
%     for j = i+1:num_groups
%         group1 = all_entropies(i);
%         group2 = all_entropies(j);
%         
%         % Perform the Mann-Whitney U test
%         [p_value, ~] = ranksum(group1, group2);
%         
%         % Check if the p-value is less than 0.05
%         if p_value < 0.05
%             significant_pairs = [significant_pairs; i, j];
%         end
%     end
% end
% % Display the significant pairs
% disp('Significant pairs (p < 0.05):');
% disp(significant_pairs);

avg_values_cached = avg_values;
figure;
hold on;
for i = 1:numel(avg_values)
    if i == 2
        avg_values(i) = avg_values(i) / 6;
    elseif i == 3
        avg_values(i) = avg_values(i) / 8;
    elseif i ==4
        avg_values(i) = avg_values(i) / 5;
    end
    bar_handle = bar(i, avg_values(i), 'FaceColor', colors{i});
end
hold off;

ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};

% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

avg_values = avg_values_cached;
figure;
hold on;
for i = 1:numel(avg_values)
    if i==1 || i ==5
        continue;
    end
    if i == 2
        avg_values(i) = (avg_values(i) + avg_values(1)) / 7;
    elseif i == 3
        avg_values(i) = avg_values(i) / 8;
    elseif i ==4
        avg_values(i) = (avg_values(i) + avg_values(5)) / 6;
    end
    bar_handle = bar(i-1, avg_values(i), 'FaceColor', colors{i});
end
hold off;

ordered_keys = {'HCPremanifest', 'Mild', 'SevereHGPS'};

% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

new_group_definitions = {
    [1, 2, 3, 4, 5, 6, 7],  ...
    [8, 9, 10, 11, 12, 13, 14, 15],   ...
    [16, 17, 18, 19, 20, 21]  
};
num_new_groups = length(new_group_definitions);
combined_entropies = cell(num_new_groups, 1);
for i = 1:num_new_groups
    original_group_indices = new_group_definitions{i};
    combined_entropies{i} = all_entropies(original_group_indices);
end
significant_pairs = [];
p_values = [];
for i = 1:num_new_groups
    for j = i+1:num_new_groups
        group1 = combined_entropies{i};
        group2 = combined_entropies{j};
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1, group2);
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            p_values = [p_values;p_value];
        end
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% Shannon per Cluster
close all
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

strNames = string(clusterNames);
values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
% dictCached = containers.Map(strNames, values);

currTypes = uniqueTempCAP(1:7);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p; 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end

% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{2});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('HC-Premanifest');
ylim([0 1.5]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(8:15);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p; 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end

% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{3});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Mild');
ylim([0 1.5]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(16:21);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p; 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end

% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{4});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Severe-HGPS');
ylim([0 1.5]);


% P Value
group_entropies = cell(1, length(values));
for i = 1:numel(values)
    group_entropies{i} = values(i);
end
all_entropies = vertcat(group_entropies{:});

% new_group_definitions = {
%     [1, 2, 3, 4, 5, 6, 7],  ...
%     [8, 9, 10, 11, 12, 13, 14, 15],   ...
%     [16, 17, 18, 19, 20, 21]  
% };
% num_new_groups = length(new_group_definitions);
% combined_entropies = cell(num_new_groups, 1);
% for i = 1:num_new_groups
%     original_group_indices = new_group_definitions{i};
%     combined_entropies{i} = all_entropies(original_group_indices);
% end
% significant_pairs = [];
% p_values = [];
% for i = 1:num_new_groups
%     for j = i+1:num_new_groups
%         group1 = combined_entropies{i};
%         group2 = combined_entropies{j};
%         
%         % Perform the Mann-Whitney U test
%         [p_value, ~] = ranksum(group1, group2);
%         
%         % Check if the p-value is less than 0.05
%         if p_value < 0.05
%             significant_pairs = [significant_pairs; i, j];
%             p_values = [p_values;p_value];
%         end
%     end
% end
% % Display the significant pairs
% disp('Significant pairs (p < 0.05):');
% disp(significant_pairs);
% disp(p_values);
%% Shannon per Cluster all
close all
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

strNames = string(clusterNames);
values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
% dictCached = containers.Map(strNames, values);
meansValues = [];
currTypes = uniqueTempCAP(1);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
group_entropies = cell(5, length(values));
for i = 1:numel(values)
    group_entropies{1,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{1});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('HC');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(2:7);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{2,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{2});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Premanifest');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(8:15);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{3,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{3});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Mild');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(16:20);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{4,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{4});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Severe');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(21);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{5,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{5});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('HGPS');
ylim([0 1]);


% Create the MEAN bar plot with different colors
figure;
hold on;
for i = 1:numel(meansValues)
    bar_handle = bar(i, meansValues(i), 'FaceColor', colors{i});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

% P Value
significant_pairs = [];
p_values = [];
for i = 1:5
    for j = i+1:5
        group1 = cell2mat(group_entropies(i,:));
        group2 = cell2mat(group_entropies(j,:));
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1, group2);
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            
        end
        p_values = [p_values;p_value];
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% Gini-Simpson per Cluster all
close all
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

strNames = string(clusterNames);
values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
% dictCached = containers.Map(strNames, values);
meansValues = [];
currTypes = uniqueTempCAP(1);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p^2; 
        end
    end
end
for key = keys(dict)
    dict(key{1}) = dict(key{1});
end

ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
group_entropies = cell(5, length(values));
for i = 1:numel(values)
    group_entropies{1,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{1});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('HC');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(2:7);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p^2; 
        end
    end
end
for key = keys(dict)
    dict(key{1}) = dict(key{1});
end

ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{2,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{2});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Premanifest');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(8:15);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p^2; 
        end
    end
end
for key = keys(dict)
    dict(key{1}) = dict(key{1});
end

ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{3,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{3});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Mild');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(16:20);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p^2; 
        end
    end
end
for key = keys(dict)
    dict(key{1}) = dict(key{1});
end

ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{4,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{4});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('Severe');
ylim([0 1]);

values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(21);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/clusterSize;
            dict(strNames(myGroup)) = dict(strNames(myGroup)) + p^2; 
        end
    end
end
for key = keys(dict)
    dict(key{1}) = dict(key{1});
end

ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{5,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{5});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title('HGPS');
ylim([0 1]);


% Create the MEAN bar plot with different colors
figure;
hold on;
for i = 1:numel(meansValues)
    bar_handle = bar(i, meansValues(i), 'FaceColor', colors{i});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

% P Value
significant_pairs = [];
p_values = [];
for i = 1:5
    for j = i+1:5
        group1 = cell2mat(group_entropies(i,:));
        group2 = cell2mat(group_entropies(j,:));
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1', group2');
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            
        end
        p_values = [p_values;p_value];
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% TEST Shannon per Cluster all
close all
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

strNames = string(clusterNames);

% dictCached = containers.Map(strNames, values);
meansValues = [];
group_entropies = cell(length(uniqueTempCAP), length(values));

for CAPNum=1:length(uniqueTempCAP)
    values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(CAPNum);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currCAP = sum(string(tempCAP2) == currTypes);
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/(clusterSize);
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(end+1) = mean(values);
for i = 1:numel(values)
    group_entropies{CAPNum,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;

for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{1});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title(currTypes);
ylim([0 1]);
end



% Create the MEAN bar plot with different colors
figure;
hold on;
avg_values = [];
for i = 1:numel(meansValues)
    currTypes = uniqueTempCAP(i);
    currCAP = sum(string(tempCAPCached) == currTypes);
    
    if i == 1
      col = colors{1};
      avg_values(end+1) = meansValues(i)/currCAP;
    elseif i < 8 % 6
        if i == 2
            avg_values(end+1) = meansValues(i)/currCAP;
        else
          avg_values(end) = avg_values(end) + meansValues(i)/currCAP;
        end
        col = colors{2};
    elseif i < 16 % 8
        if i == 8
            avg_values(end+1) = meansValues(i)/currCAP;
        else
          avg_values(end) = avg_values(end) + meansValues(i)/currCAP;
        end
        col = colors{3};
    elseif i < 21 % 5
        if i == 16
            avg_values(end+1) = meansValues(i)/currCAP;
        else
          avg_values(end) = avg_values(end) + meansValues(i)/currCAP;
        end
        col = colors{4};
    elseif i == 21
        avg_values(end+1) = meansValues(i)/currCAP;
        col = colors{5};
    end
    bar_handle = bar(i, meansValues(i)/currCAP, 'FaceColor', col);
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
ordered_keys = uniqueTempCAP;
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

figure;
hold on;
for i = 1:numel(avg_values)
    if i == 2
        avg_values(i) = avg_values(i) / 6;
    elseif i == 3
        avg_values(i) = avg_values(i) / 8;
    elseif i ==4
        avg_values(i) = avg_values(i) / 5;
    end
    bar_handle = bar(i, avg_values(i), 'FaceColor', colors{i});
end
hold off;

ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};

% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);

% P Value
significant_pairs = [];
p_values = [];
for i = 1:length(uniqueTempCAP)
    for j = i+1:length(uniqueTempCAP)
        group1 = cell2mat(group_entropies(i,:));
        group2 = cell2mat(group_entropies(j,:));
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1, group2);
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            p_values = [p_values;p_value];
        end
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% TEST2 Shannon per Cluster all
close all
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

strNames = string(clusterNames);
strNames = strNames(1:end-1);
% dictCached = containers.Map(strNames, values);
meansValues = cell(21, 6);
group_entropies = cell(length(uniqueTempCAP), numel(strNames));
for CAPNum=1:length(uniqueTempCAP)
    values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(CAPNum);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currCAP = sum(string(tempCAPCached) == currTypes);
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/(clusterSize);
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p))/currCAP;
%             p = p / currCAP;
%             dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(CAPNum,:) = num2cell(values);
for i = 1:numel(values)
    group_entropies{CAPNum,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;

for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{1});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title(currTypes);
% ylim([0 16e-4]);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, strcat('\Shannon_',string(CAPNum),'.png')))
end

meansValuesMat= cell2mat(meansValues);
groups = {1, 2:7, 8:15, 16:20, 21};
groupMeansMat = zeros(numel(groups), size(meansValuesMat, 2));
for i = 1:numel(groups)
    group = groups{i};
    if numel(group) == 1
        groupMeansMat(i, :) = meansValuesMat(group, :);
    else
        groupMeansMat(i, :) = mean(meansValuesMat(group, :));
    end
end

figure();
meansMat = mean(meansValuesMat, 2);
for i = 1:numel(meansMat)
    if i == 1
      col = colors{1};
    elseif i < 8 % 6
        col = colors{2};
    elseif i < 16 % 8
        col = colors{3};
    elseif i < 21 % 5
        col = colors{4};
    elseif i == 21
        col = colors{5};
    end
    bar_handle = bar(i, meansMat(i), 'FaceColor', col);
    hold on;
end
p = polyfit(1:numel(meansMat), meansMat, 3);
y_fit = polyval(p, 1:numel(meansMat));
plot(1:numel(meansMat), y_fit, 'k', 'LineWidth', 2);
y_mean = mean(meansMat);
SSR = sum((y_fit - y_mean).^2);
SST = sum((meansMat - y_mean).^2);
R2 = SSR / SST;
% text(numel(meansMat) * 0.6, max(meansMat) * 0.9, sprintf('R^2 = %.2f', R2), 'FontSize', 12);
hold off;

% Set the Y-axis label and X-axis tick labels
ordered_keys2 = uniqueTempCAP;
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys2));
ordered_keys2(1) = 'HC';
ordered_keys2(end) = 'HGPS';
xticklabels(ordered_keys2);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, strcat('\Shannon_All_1.png')))
% Create the MEAN bar plot with different colors
% Iterate through the rows of groupMeansMat
for i = 1:size(groupMeansMat, 1)
    figure;
    bar(groupMeansMat(i, :),'FaceColor', colors{i});  % Plot a bar graph for the current row
    set(gca, 'XTickLabel', ordered_keys, 'XTickLabelRotation', 45);  % Set x-axis labels and rotate them by 45 degrees
    if i == 1
        title('HC');
    elseif i == 2
        title('Premanifest');
    elseif i == 3
        title('Mild');
    elseif i == 4
        title('Severe'); 
    elseif i == 5
        title('HGPS');
    end
            
    ylabel('Shannon Entropy');
    ylim([0 16e-4]);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, strcat('\Shannon_',string(i),'_b.png')))
end

figure();
rowMeans = mean(groupMeansMat, 2);
hold on
for i = 1:numel(rowMeans)
    bar_handle = bar(i, rowMeans(i), 'FaceColor', colors{i});
end
p = polyfit(1:numel(rowMeans), rowMeans, 3);
y_fit = polyval(p, 1:numel(rowMeans));
plot(1:numel(rowMeans), y_fit, 'k', 'LineWidth', 2);
y_mean = mean(rowMeans);
SSR = sum((y_fit - y_mean).^2);
SST = sum((rowMeans - y_mean).^2);
R2 = SSR / SST;
% text(numel(rowMeans) * 0.6, max(rowMeans) * 0.9, sprintf('R^2 = %.2f', R2), 'FontSize', 12);
hold off;

% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
% ylim([0 16e-4]);
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, strcat('\Shannon_All_2.png')))
% P Value
significant_pairs = [];
p_values = [];
p_values_all = [];

for i = 1:size(groupMeansMat,1)
    for j = i+1:size(groupMeansMat,1)
        group1 = groupMeansMat(i,:);
        group2 = groupMeansMat(j,:);
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1, group2);
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            p_values = [p_values;p_value];
        end
        p_values_all = [p_values_all;p_value];
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% TEST2 Shannon by Clusters Clusters
close all
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

% subCell = [1, 2, 6];
% subCell = [3, 4, 5];
% subCell = [7];
subCell = [1,2,3,4,5,6,7];
clusterNamesUse = clusterNames(subCell);
clusterLabelsUse = clusterLabels(subCell);
strNames = string(clusterNamesUse);

% dictCached = containers.Map(strNames, values);
meansValues = cell(21, length(subCell));
group_entropies = cell(length(uniqueTempCAP), numel(strNames));

for CAPNum=1:length(uniqueTempCAP)
    values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(CAPNum);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currCAP = sum(string(tempCAPCached) == currTypes);
        currentLookedGroup = groupsCell(clusterLabelsUse{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/(clusterSize);
            dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p))/currCAP;
%             p = p / currCAP;
%             dict(strNames(myGroup)) = dict(strNames(myGroup)) - (p * log(p)); 

        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(CAPNum,:) = num2cell(values);
for i = 1:numel(values)
    group_entropies{CAPNum,i} = values(i);
end
% Create the bar plot with different colors
% figure;
% hold on;

% for i = 1:numel(values)
%     bar_handle = bar(i, values(i), 'FaceColor', colors{1});
% end
% hold off;
% Set the Y-axis label and X-axis tick labels
% ylabel('Shannon Entropy');
% xticks(1:numel(ordered_keys));
% xticklabels(ordered_keys);
% title(currTypes);
% ylim([0 16e-4]);
% saveas(gcf,strcat('C:\Cell_Migration_MATLAB\NatureAviv\plts\', string(CAPNum),'.jpg'))
end

meansValuesMat= cell2mat(meansValues);
groups = {1, 2:7, 8:15, 16:20, 21};
groupMeansMat = zeros(numel(groups), size(meansValuesMat, 2));
for i = 1:numel(groups)
    group = groups{i};
    if numel(group) == 1
        groupMeansMat(i, :) = meansValuesMat(group, :);
    else
        groupMeansMat(i, :) = mean(meansValuesMat(group, :));
    end
end

figure();
meansMat = mean(meansValuesMat, 2);
for i = 1:numel(meansMat)
    if i == 1
      col = colors{1};
    elseif i < 8 % 6
        col = colors{2};
    elseif i < 16 % 8
        col = colors{3};
    elseif i < 21 % 5
        col = colors{4};
    elseif i == 21
        col = colors{5};
    end
    bar_handle = bar(i, meansMat(i), 'FaceColor', col);
    hold on;
end
p = polyfit(1:numel(meansMat), meansMat, 3);
y_fit = polyval(p, 1:numel(meansMat));
plot(1:numel(meansMat), y_fit, 'k', 'LineWidth', 2);
y_mean = mean(meansMat);
SSR = sum((y_fit - y_mean).^2);
SST = sum((meansMat - y_mean).^2);
R2 = SSR / SST;
text(numel(meansMat) * 0.6, max(meansMat) * 0.9, sprintf('R^2 = %.2f', R2), 'FontSize', 12);
hold off;

% Set the Y-axis label and X-axis tick labels
ordered_keys2 = uniqueTempCAP;
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys2));
ordered_keys2(1) = 'HC';
ordered_keys2(end) = 'HGPS';
xticklabels(ordered_keys2);
% title(strNames);
% ylim([0 36e-4]);
ylim([0 9e-4]);

% Create the MEAN bar plot with different colors
% Iterate through the rows of groupMeansMat
% for i = 1:size(groupMeansMat, 1)
%     figure;
%     bar(groupMeansMat(i, :),'FaceColor', colors{i});  % Plot a bar graph for the current row
%     set(gca, 'XTickLabel', ordered_keys, 'XTickLabelRotation', 45);  % Set x-axis labels and rotate them by 45 degrees
%     if i == 1
%         title('HC');
%     elseif i == 2
%         title('Premanifest');
%     elseif i == 3
%         title('Mild');
%     elseif i == 4
%         title('Severe'); 
%     elseif i == 5
%         title('HGPS');
%     end
%             
%     ylabel('Shannon Entropy');
%     ylim([0 16e-4]);
% 
% end

figure();
rowMeans = mean(groupMeansMat, 2);
hold on
for i = 1:numel(rowMeans)
    bar_handle = bar(i, rowMeans(i), 'FaceColor', colors{i});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
% ylim([0 36e-4]);
ylim([0 7e-4]);


% title(strNames);

% P Value
significant_pairs = [];
p_values = [];
p_values_all = [];

for i = 1:size(groupMeansMat,1)
    for j = i+1:size(groupMeansMat,1)
        group1 = groupMeansMat(i,:);
        group2 = groupMeansMat(j,:);
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1, group2);
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            p_values = [p_values;p_value];
        end
        p_values_all = [p_values_all;p_value];
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% TEST2 Gini-Simpson per Cluster all
close all
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
colors = newcolorsRGB; % [blue; green; orange; red; purple]

tempCAP = data.CAP(combinedIndex);
tempCAP(string(tempCAP) == "76.42527") = "77";
uniqueTempCAP = string(unique(tempCAP));
numericTempCAP = str2double(uniqueTempCAP);
roundedNumericTempCAP = round(numericTempCAP);
uniqueTempCAP = string(roundedNumericTempCAP);

tempCAPCached = tempCAP;
roundedNumericTempCAPCached = round(tempCAPCached);
tempCAPCached = string(roundedNumericTempCAPCached);

strNames = string(clusterNames);

% dictCached = containers.Map(strNames, values);
meansValues = cell(21, 7);
group_entropies = cell(length(uniqueTempCAP), length(values));

for CAPNum=1:length(uniqueTempCAP)
    values = zeros(1, numel(strNames));
dict = containers.Map(strNames, values);
currTypes = uniqueTempCAP(CAPNum);
for type=1:length(currTypes)
    for myGroup=1:1:length(strNames)
        currCAP = sum(string(tempCAPCached) == currTypes);
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempCAP(logical(aNumbers));
        numericTempCAP = str2double(string(types));
        roundedNumericTempCAP = round(numericTempCAP);
        types = string(roundedNumericTempCAP);
        clusterSize = length(types);
        num = sum(types == currTypes(type));
        if num > 0
           p = num/(clusterSize);
            dict(strNames(myGroup)) = dict(strNames(myGroup)) +(p^2)/currCAP; 
        end
    end
end


ordered_keys = strNames;

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
valuesCached = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict(key);
end
meansValues(CAPNum,:) = num2cell(values);
for i = 1:numel(values)
    group_entropies{CAPNum,i} = values(i);
end
% Create the bar plot with different colors
figure;
hold on;

for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{1});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
title(currTypes);
ylim([0 0.00006]);
end

meansValuesMat= cell2mat(meansValues);
groups = {1, 2:7, 8:15, 16:20, 21};
groupMeansMat = zeros(numel(groups), size(meansValuesMat, 2));
for i = 1:numel(groups)
    group = groups{i};
    if numel(group) == 1
        groupMeansMat(i, :) = meansValuesMat(group, :);
    else
        groupMeansMat(i, :) = mean(meansValuesMat(group, :));
    end
end

% Create the MEAN bar plot with different colors
% Iterate through the rows of groupMeansMat
for i = 1:size(groupMeansMat, 1)
    figure;
    bar(groupMeansMat(i, :),'FaceColor', colors{i});  % Plot a bar graph for the current row
    set(gca, 'XTickLabel', ordered_keys, 'XTickLabelRotation', 45);  % Set x-axis labels and rotate them by 45 degrees
    if i == 1
        title('HC');
    elseif i == 2
        title('Premanifest');
    elseif i == 3
        title('Mild');
    elseif i == 4
        title('Severe'); 
    elseif i == 5
        title('HGPS');
    end
            
    ylabel('Gini-Simpson Index');
ylim([0 0.00006]);

end

figure();
rowMeans = mean(groupMeansMat, 2);
hold on
for i = 1:numel(rowMeans)
    bar_handle = bar(i, rowMeans(i), 'FaceColor', colors{i});
end
hold off;
% Set the Y-axis label and X-axis tick labels
ylabel('Gini-Simpson Index');
ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
ylim([0 0.00006]);

% P Value
significant_pairs = [];
p_values = [];
p_values_all = [];

for i = 1:size(groupMeansMat,1)
    for j = i+1:size(groupMeansMat,1)
        group1 = groupMeansMat(i,:);
        group2 = groupMeansMat(j,:);
        
        % Perform the Mann-Whitney U test
        [p_value, ~] = ranksum(group1, group2);
        
        % Check if the p-value is less than 0.05
        if p_value < 0.05
            significant_pairs = [significant_pairs; i, j];
            p_values = [p_values;p_value];
        end
        p_values_all = [p_values_all;p_value];
    end
end
% Display the significant pairs
disp('Significant pairs (p < 0.05):');
disp(significant_pairs);
disp(p_values);
%% Shannon
tempType = data.CellType(combinedIndex);
tempTypeCached = tempType;
uniqueTempType = string(unique(tempType));
dict = struct();
% dict.HCPremanifest = 0;
% dict.Mild = 0;
% dict.SevereHGPS = 0;
dict.HC = 0;
dict.HGPS = 0;
dict.Mild = 0;
dict.Premanifest = 0;
dict.Severe = 0;

counter = 0;
for myGroup=1:1:length(clusterLabels)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempType(logical(aNumbers));
% %         types(string(types) == "'HC'") = "HCPremanifest";
%         types(string(types) == "'Premanifest'") = "HCPremanifest";
%         types(string(types) == "'Severe'")    = "SevereHGPS";
%         types(string(types) == "'HGPS'") = "SevereHGPS";

        num_all = length(types);
        uniqueTypes = unique(types);
        for type=1:length(uniqueTypes)
            num_all_new = sum(string(tempType) == string(uniqueTypes(type)));
            group_size = length(tempTypeCached(tempTypeCached == uniqueTypes(type)));
            num = sum(string(types) == string(uniqueTypes(type)));
            p = num/num_all_new;
%             if type == 3
%             disp(num_all);
%             disp(num);
%             counter = counter + num;
%             end
            
            name = string(uniqueTypes(type));
            name = name{1};
%             if name == "'Mild'"
                name = name(2:end-1);
%             end
%             if name == "HC" || name == "Premanifest"
%                 name = "HCPremanifest";
%             elseif name == "Severe" || name == "HGPS"
%                 name = "SevereHGPS";
%             end
            
            dict.(name) = dict.(name) - (p * log(p));
%             dict.(name) = dict.(name) + p;

        end
end
% Define the order and colors
% ordered_keys = {'HCPremanifest', 'Mild', 'SevereHGPS'};
ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};

colors = newcolorsRGB; % [blue; green; orange; red; purple]

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict.(key);
end

% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{i});
end
hold off;

% Set the Y-axis label and X-axis tick labels
ylabel('Shannon Entropy');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
%% Simpson
tempType = data.CellType(combinedIndex);
uniqueTempType = string(unique(tempType));
dict = struct();
% dict.HCPremanifest = 0;
% dict.Mild = 0;
% dict.SevereHGPS = 0;
dict.HC = 0;
dict.HGPS = 0;
dict.Mild = 0;
dict.Premanifest = 0;
dict.Severe = 0;

for myGroup=1:1:length(clusterLabels)
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        currNums = str2double(regexp(currentLookedGroup,'\d*','Match', 'once'));
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(currNums) = 1;
        types = tempType(logical(aNumbers));
%         types(string(types) == "'HC'") = "HCPremanifest";
%         types(string(types) == "'Premanifest'") = "HCPremanifest";
%         types(string(types) == "'Severe'") = "SevereHGPS";
%         types(string(types) == "'HGPS'") = "SevereHGPS";

        num_all = length(types);
        uniqueTypes = unique(types);
        for type=1:length(uniqueTypes)
            num_all_new = sum(string(tempType) == string(uniqueTypes(type)));
            num = sum(string(types) == string(uniqueTypes(type)));
            p = num/num_all_new;
            name = string(uniqueTypes(type));
            name = name{1};
%             if name == "'Mild'"
                name = name(2:end-1);
%             end
%             if name == "HC" || name == "Premanifest"
%                 name = "HCPremanifest";
%             elseif name == "Severe" || name == "HGPS"
%                 name = "SevereHGPS";
%             end
            dict.(name) = dict.(name) + p^2;
        end
end

% Calculate Simpson index
% for type = 1:length(uniqueTypes)
%     name = string(uniqueTypes(type));
%     name = name{1};
%     if name == "'Mild'"
%         name = name(2:end-1);
%     end
% 
%     dict.(name) = 1 - dict.(name);
% end

% Define the order and colors
% ordered_keys = {'HCPremanifest', 'Mild', 'SevereHGPS'};
ordered_keys = {'HC', 'Premanifest', 'Mild', 'Severe', 'HGPS'};

colors = newcolorsRGB; % [blue; green; orange; red; purple]

% Extract values in the specified order
values = zeros(1, numel(ordered_keys));
for i = 1:numel(ordered_keys)
    key = ordered_keys{i};
    values(i) = dict.(key);
end

% Create the bar plot with different colors
figure;
hold on;
for i = 1:numel(values)
    bar_handle = bar(i, values(i), 'FaceColor', colors{i});
end
hold off;

% Set the Y-axis label and X-axis tick labels
ylabel('Simpson Index');
xticks(1:numel(ordered_keys));
xticklabels(ordered_keys);
%%
% Prepare 1-dimensional displacement trajectories

xys = xysMain(idxComb);

one_dim_xys = zeros([length(xys),length(xys{1})-1]);
one_dim_xys_not_normalized = zeros([length(xys),length(xys{1})-1]);

for i=1:length(xys)
    temp_traj = xys{i};
    temp_one_dim_traj = sqrt(diff(temp_traj(:,1),1).^2 + diff(temp_traj(:,2),1).^2);
%     normalized = (temp_one_dim_traj - mean(temp_one_dim_traj))./std(temp_one_dim_traj);
    one_dim_xys_not_normalized(i,:) = temp_one_dim_traj';
    normalized = zscore(temp_one_dim_traj);
    one_dim_xys(i,:) = normalized';
end
% one_dim_xys2 = zscore(one_dim_xys_not_normalized, 2);

r = randi([0 length(xys)],1,1);
xy = xys{r};
xy1 = smooth(xy(:,1));
xy2 = smooth(xy(:,2));
sz = 50;
c = linspace(1,10,length(xy1));
figure()
scatter(xy1,xy2,sz, c, 'filled'); hold on;
colorbar
colormap jet
plot(xy1,xy2,'-','color','k', 'linewidth',0.5); hold off;
title('X-Y Trajectory')
xlabel('X') 
ylabel('Y') 
figure()
xAxis = 0:10:10*73;
c = linspace(1,10,length(xAxis));
scatter(xAxis, one_dim_xys(r,:),sz, c, 'filled'); hold on;
colorbar
colormap jet
plot(xAxis, one_dim_xys(r,:),'-','color','k', 'linewidth',0.5); hold off;
title('Activity Profile')
xlabel('Elapsed Time (750 minutes)') 
ylabel('Z-Score Displacements') 


%
x_data = one_dim_xys;

% Create clustering algorithm function
myfunc = @(X,K)clusterdata(X,'Distance','cityblock','MaxClust',K, 'linkage',...
    'ward');

% evaluate optimal number of clusters
eva = evalclusters(x_data,myfunc,'CalinskiHarabasz','KList',...
    [1:40]);

display(strcat('Optimal number of clusters is', {' '} ,num2str(eva.OptimalK)));


maxclust = eva.OptimalK;

Cellular_Activity_clusters = clusterdata(x_data,'Distance','cityblock','MaxClust',maxclust,'linkage',...
    'ward');
%% Activity for Groups
close all
numOfRandoms = 10;
sz = 10;

temp_cellType = DataFileWell.CellType(idxComb);

a1Pol_Activity = x_data(temp_cellType == "'HC'",:);
a2Pol_Activity = x_data(temp_cellType == "'Severe'",:);
a3Pol_Activity = x_data(temp_cellType == "'Mild'",:);
a4Pol_Activity = x_data(temp_cellType == "'Premanifest'",:);
a5Pol_Activity = x_data(temp_cellType == "'HGPS'",:);

classes = {a1Pol_Activity, a4Pol_Activity, a3Pol_Activity, a2Pol_Activity, a5Pol_Activity};
classNames = {'HC', 'Premanifest', 'Mild', 'Severe',  'HGPS'};
for i = 1:length(classes)
    aPol = classes{i};
    currNums = [];
        figure()
        mover = 0;
        for j=1:numOfRandoms
            currNum = round(size(aPol,1)*rand());
            while ismember(currNum, currNums) || currNum == 0
                currNum = ceil(size(aPol,1)*rand());
            end
            xy=aPol(currNum,:); 
            xy = xy';
            xAxis = 0:10:10*73;
            c = linspace(1,10,length(xAxis));
            scatter(xAxis, xy + mover,sz, c, 'filled'); hold on;
            colorbar
            colormap jet
            plot(xAxis, xy + mover,'-','color','k', 'linewidth',0.5); hold on;
            title(['Activity Profile for ' classNames{i}])
            xlabel('Elapsed Time (750 minutes)') 
            ylabel('Z-Score Displacements') 
            set(gca, 'YTick', [])
            mover = mover + 10;
        end
    g = gcf;
    exportgraphics(g, strcat(folderCurrent, '\Ag_Activity_Group_' + string(i) + '.png'), 'Resolution', 900)
end
%% Lags and Trains for Groups
% binarize
% close all
meanArray = mean(one_dim_xys,2);
stdArray = std(one_dim_xys,0,2);
binarized = one_dim_xys>meanArray+0.25*stdArray;
TrainSum = [];
LagsSum = [];
TrainCV = [];
LagsCV = [];
grp1 = [];
labels = [];


a1Pol_Activity = binarized(temp_cellType == "'HC'",:);
a2Pol_Activity = binarized(temp_cellType == "'Severe'",:);
a3Pol_Activity = binarized(temp_cellType == "'Mild'",:);
a4Pol_Activity = binarized(temp_cellType == "'Premanifest'",:);
a5Pol_Activity = binarized(temp_cellType == "'HGPS'",:);

classes = {a1Pol_Activity, a4Pol_Activity, a3Pol_Activity, a2Pol_Activity, a5Pol_Activity};
classNames = {'HC', 'Premanifest', 'Mild', 'Severe',  'HGPS'};

for i = 1:length(classes)
    aPol = classes{i};
    TrainsPol = aPol;
    LagsPol = 1-aPol;
    TrainSum = [TrainSum, sum(TrainsPol, 2)'];
    LagsSum = [LagsSum, sum(LagsPol, 2)'];
    TrainCV = [TrainCV, (std(TrainsPol,0, 2)./mean(TrainsPol, 2))'];
    LagsCV = [LagsCV, (std(LagsPol,0, 2)./mean(LagsPol, 2))'];
    grp1 = [grp1, i*ones(1,length(TrainsPol))];
    
end

figure()
h = boxplot(TrainSum, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', classNames);
title('Boxplot for Trains');
ylabel('Train Length');
ylim([0 max(TrainSum)*1.005]);
medians = grpstats(TrainSum, grp1, 'mean');
x_locations = 1:numel(medians);
hold on;
plot(x_locations, medians, '-m');
hold off;
g = gcf;
exportgraphics(g, strcat(folderCurrent, '\Ah_Groups_Trains_Length.png'), 'Resolution', 900)
figure()
boxplot(LagsSum, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', classNames)
title('Boxplot for Lags');
ylabel('Lags Length');
ylim([0 max(LagsSum)*1.05]);
medians = grpstats(LagsSum, grp1, 'mean');
x_locations = 1:numel(medians);
hold on;
plot(x_locations, medians, '-m');
hold off;
g = gcf;
exportgraphics(g, strcat(folderCurrent, '\Ah_Groups_Lags_Length.png'), 'Resolution', 900)
figure()
boxplot(TrainCV, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', classNames)
title('Boxplot for CV Trains');
ylabel('CV Train Length');
ylim([0 max(TrainCV)*1.05]);
medians = grpstats(TrainCV, grp1, 'mean');
x_locations = 1:numel(medians);
hold on;
plot(x_locations, medians, '-m');
hold off;
g = gcf;
exportgraphics(g, strcat(folderCurrent, '\Ah_Groups_CVTrains_Length.png'), 'Resolution', 900)
figure()
boxplot(LagsCV, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', classNames)
title('Boxplot for CV Lags');
ylabel('CV Lags Length');
ylim([0 max(LagsCV)*1.05]);
medians = grpstats(LagsCV, grp1, 'mean');
x_locations = 1:numel(medians);
hold on;
plot(x_locations, medians, '-m');
hold off;
g = gcf;
exportgraphics(g, strcat(folderCurrent, '\Ah_Groups_CVLags_Length.png'), 'Resolution', 900)

%% Clusters
close all force
runInsiders = false;
if runInsiders
    folderCurrent = folderCurrent + 'Insiders\';
    mkdir(folderCurrent);
end
cMap = [linspace(0,1,128)'*[1,1], ones(128,1)]; % The blue range
cMap = [cMap; rot90(cMap,2)];                   % add the red range
cMap = (cMap);

temp_cellType = DataFileWell.CellType(idxComb);
temp_CAPPatient = string(DataFileWell.CAPPatient(idxComb));

rowLabels = string(temp_cellType);
for i=1:size(rowLabels, 1)
    rowLabels(i) = rowLabels(i) + '_' + string(i) + '_' + temp_CAPPatient(i);
end
cg = clustergram(one_dim_xys, ...
                             'RowLabels', rowLabels,...
                             'RowPdist', 'cityblock',...
                             'ColumnPdist', 'cityblock',...
                             'Linkage', 'ward',...
                             'Cluster', 'column',...
                             'Dendrogram', 0,'Colormap', cMap);
%
% GI_One_Dim = [go1; go2; go3; go4; go5; go6; go7];
% GI_One_Dim = [g1; g2; g5];
% GI_One_Dim = [g1; g2; g3; g4];
GINums = [6003 6000 5996 6002];
% GINums = [6003 6000 6004];
% GINums = [5982 5995 5987 5994 5999 6001 5996]; good
% cluster Group
de = 1.1;
seperate = 0; %33 32 42 41 40

clear groupsCell
numGroups = length(cg.RowLabels) - 1;
numGroups = length(GI_One_Dim);
clusterDE = {};
clusterHC = {};
clusterR = {};
clusterM = {};
clusterS = {};
clusterP = {};
clusterLabels = {};
clusterNames = {};
myCellHC = {};
myCellR = {};
myCellM = {};
myCellS = {};
myCellP = {};
myCellHCCAPS = {};
myCellRCAPS = {};
myCellMCAPS = {};
myCellSCAPS = {};
myCellPCAPS = {};

% groupsCell = cell([numGroups 1]);
groupIndexed=0;
uniqueCAPS = natsort(unique(regexp(string(data.CAPPatient2(combinedIndex)),'\d*','Match', 'once')));
for groupI=numGroups:-1:1
    group = GINums(groupI);
    myBreakFlag = false;
    flag = false;
    currentGroupInfo = getGroupInfo(cg,group, 1);
    currentGroupInfo = GI_One_Dim(groupI);
    groupsCell(group) = currentGroupInfo;
    if ~runInsiders
        for iterGroup=1:length(clusterLabels)
            currNames = groupsCell(clusterLabels{iterGroup}).RowNodeNames;
            if sum(ismember(currNames, currentGroupInfo.RowNodeNames)) ~= 0
                myBreakFlag = true;
                break;
            end
        end
        if myBreakFlag
           continue; 
        end
    end
    countSevere = sum(contains(string(currentGroupInfo.RowNodeNames), "Severe"));
    countMild = sum(contains(string(currentGroupInfo.RowNodeNames), "Mild"));
    countRisk = sum(contains(string(currentGroupInfo.RowNodeNames), "Premanifest"));
    countHC = sum(contains(string(currentGroupInfo.RowNodeNames), "HC"));
    countProgeria = sum(contains(string(currentGroupInfo.RowNodeNames), "HGPS"));
    
    CAPScores = containers.Map;
    for uC=1:length(uniqueCAPS)
        CAPScores(uniqueCAPS(uC)) = 0;
    end
    for item=1:length(currentGroupInfo.RowNodeNames)
        currString = string(currentGroupInfo.RowNodeNames(item));
        split_string = strsplit(currString, '_');
        resultCAP = regexp(split_string{3},'\d*','Match', 'once');
        CAPScores(resultCAP) = CAPScores(resultCAP) + 1;
    end
    [~,sortedCAPScoresIndexes] = natsort(CAPScores.keys);
    CAPScoresString = string(CAPScores.values);
    CAPScoresString = CAPScoresString(sortedCAPScoresIndexes);
    
    if countMild < countHC && countRisk < countHC && countSevere < countHC && countProgeria < countHC% HC
    if max([countMild countRisk countSevere countProgeria]) / countHC < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countMild countRisk countSevere countProgeria]) / countHC};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    if countMild < countSevere && countRisk < countSevere && countHC < countSevere && countProgeria < countSevere% Severe
    if max([countMild countRisk countHC countProgeria]) / countSevere < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countMild countRisk countHC countProgeria]) / countSevere};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end

    if countSevere < countMild && countRisk < countMild && countHC < countMild && countProgeria < countMild% Mild
    if max([countSevere countRisk countHC countProgeria]) / countMild < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countSevere countRisk countHC countProgeria]) / countMild};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    
    if countSevere < countRisk && countMild < countRisk && countHC < countRisk && countProgeria < countRisk% Risk
    if max([countSevere countMild countHC countProgeria]) / countRisk < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countSevere countMild countHC countProgeria]) / countRisk};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    
    if countSevere < countProgeria && countMild < countProgeria && countHC < countProgeria && countRisk < countProgeria% Progeria
    if max([countSevere countMild countHC countRisk]) / countProgeria < de
        if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
            flag = true;
            numbers = str2double(regexp(string(currentGroupInfo.RowNodeNames),'\d*','Match', 'once'));
            if sum(ismember(numbers, group)) ~= 0
%                continue; 
            end
            clusterLabels(end+1) = {group};
            clusterDE(end+1) = {max([countSevere countMild countHC countRisk]) / countProgeria};
            clusterHC(end+1) = {countHC};
            clusterR(end+1) = {countRisk};
            clusterM(end+1) = {countMild};
            clusterS(end+1) = {countSevere};
            clusterP(end+1) = {countProgeria};
        end
    end
    end
    
        if flag
            disp("******************************************************");
            checkFlag = false;
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countHC
            checkFlag = true;
           disp(strcat(string(group), "_HC_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed))};
                myCellHC = [myCellHC; strcat('Cluster-',string(groupIndexed)), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellHCCAPS = [myCellHCCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];

           end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countSevere
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_Severe_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed))};
                myCellS = [myCellS; strcat('Cluster-',string(groupIndexed)), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellSCAPS = [myCellSCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];
           end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countMild
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_Mild_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_Progeria_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed))};
                myCellM = [myCellM; strcat('Cluster-',string(groupIndexed)), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellMCAPS = [myCellMCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];

           end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countRisk
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_Risk_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed))};
                myCellR = [myCellR; strcat('Cluster-',string(groupIndexed)), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellRCAPS = [myCellRCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];
            end
        end
        if max([countSevere, countMild, countRisk, countHC, countProgeria]) ==  countProgeria
            if checkFlag
                 disp("ERROR&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
            end
            checkFlag = true;
           disp(strcat(string(group), "_HGPS_Cluster"));
           disp(strcat("Num_Of_HC_", string(countHC)));
           disp(strcat("Num_Of_Severe_", string(countSevere)));
           disp(strcat("Num_Of_Mild_", string(countMild)));
           disp(strcat("Num_Of_Premanifest_", string(countRisk)));
           disp(strcat("Num_Of_HGPS_", string(countProgeria)));
           if (countMild + countRisk + countSevere + countHC + countProgeria) > seperate
               groupIndexed = groupIndexed + 1;
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed))};
                myCellP = [myCellP; strcat('Cluster-',string(groupIndexed)), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellPCAPS = [myCellPCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];
            end
        end
        
        disp(strcat("Total_Cells_Is_", string(countMild + countRisk + countSevere + countHC + countProgeria)));

        end
end
    
groupsCellActivity = groupsCell;
clusterNamesActivity = clusterNames;
%
CM = jet(length(clusterLabels));
CM(:,4) = 0.5;
CMCell = {};
for color=1:size(CM,1)
    CMCell(end+1) = {CM(color,:)};
end
CMCell = CMCell(randperm(numel(CMCell)));

%
if ~isempty(clusterLabels)
    cNames = clusterNames;
    clusterNames{1} = "Cluster-3";
    clusterNames{2} = "Cluster-2";
    clusterNames{3} = "Cluster-4";
    clusterNames{4} = "Cluster-1";
rm = struct('GroupNumber',clusterLabels,'Annotation',clusterNames,...
    'Color',CMCell, 'FontSize', 5);
set(cg,'RowGroupMarker',rm)
cgf = plot(cg);
cgf = plot(cg);
set(cgf,'FontSize',14, 'fontweight','bold')
fig = get(cgf, 'Parent');
fig.WindowState = 'maximized';
saveas(cgf,strcat(folderCurrent, '\c_Activity_Clustergram.jpg'))
end
%%
% close all force

%$%$ BAR PLOT
myCellHC(1) = "Cluster-1";
myCellM(1) = "Cluster-2";
myCellR(1) = "Cluster-3";
myCellP(1) = "Cluster-4";
myCellTotal = [myCellHC; myCellM; myCellR; myCellS; myCellP];
XAxis = categorical(myCellTotal(:,1));
XAxis = reordercats(XAxis,flip(myCellTotal(:,1)));
YAxis = str2double(myCellTotal(:,2:end));
sumY = sum(YAxis, 2);
Answer = bsxfun(@rdivide,YAxis,sumY(:));

YAxis1 = Answer(:,1);
YAxis2 = Answer(:,2);
YAxis3 = Answer(:,3);
YAxis4 = Answer(:,4);
YAxis5 = Answer(:,5);

Xmax = 1.05;

figure();
barh(XAxis, YAxis1);
colororder(newcolors{1})
title('HC Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\Ad_Bar_1.png'))
figure();
barh(XAxis, YAxis2);
colororder(newcolors{2})
title('Premanifest Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\Ad_Bar_2.png'))
figure();
barh(XAxis, YAxis3);
colororder(newcolors{3})
title('Mild Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\Ad_Bar_3.png'))
figure();
barh(XAxis, YAxis4);
colororder(newcolors{4})
title('Severe Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\Ad_Bar_4.png'))
figure();
barh(XAxis, YAxis5);
colororder(newcolors{5})
title('HGPS Clusters');
xlim([0 Xmax])
set(gca,'FontSize',14, 'fontweight','bold')
g = gcf;
saveas(g, strcat(folderCurrent, '\Ad_Bar_5.png'))
%%
%$%$ BAR PLOT CAP
myCellTotal = [myCellHCCAPS; myCellRCAPS; myCellMCAPS; myCellSCAPS; myCellPCAPS];
XAxis = categorical(myCellTotal(:,1));
XAxis = reordercats(XAxis,flip(myCellTotal(:,1)));
YAxis = str2double(myCellTotal(:,2:end));
% sumY = sum(YAxis, 2);
% Answer = bsxfun(@rdivide,YAxis,sumY(:));

Xmax = max(0.5);
for item=1:size(myCellTotal,2)-1
%     YAxis = Answer(:,item);
    
    YAxisPlot = YAxis(:,item)/sum(YAxis(:,item));
%     Xmax = max(YAxisPlot);
    figure();
    barh(XAxis, YAxisPlot);
    colororder(newcolors{6})
    title(uniqueCAPS(item));
    if item == 1
        title('HC');
    end
    if item == size(myCellTotal,2)-1
        title('HGPS');
    end
    xlim([0 Xmax])
    set(gca,'FontSize',14, 'fontweight','bold')
    g = gcf;
    saveas(g, strcat(folderCurrent, '\Ae_Bar_' + string(item) + '.png'))
end

%% Activity
close all
numOfRandoms = 10;
sz = 10;
if ~isempty(clusterLabels)
    
    for iBig=1:1:length(clusterLabels)
    starts = iBig;
    ends = iBig+0;
    if ends > length(clusterLabels)
        ends = length(clusterLabels); 
    end
%     fig = figure('Position', get(0, 'Screensize'));
%     t = tiledlayout(3,1, 'TileSpacing','Compact');
    for myGroup=starts:ends
        currentGroupName = clusterNames{myGroup};
%         if strfind(currentGroupName, 'Premanifest') > 0
%             myColorO = newcolors{2};
%             myTypeOuter = 'Risk';
%         elseif strfind(currentGroupName, 'Mild') > 0
%             myColorO = newcolors{3};
%             myTypeOuter = 'Mild';
%         elseif strfind(currentGroupName, 'Severe') > 0
%             myColorO = newcolors{4};
%             myTypeOuter = 'Severe';
%         elseif strfind(currentGroupName, 'HC') > 0
%             myColorO = newcolors{1};
%             myTypeOuter = 'HC';
%         elseif strfind(currentGroupName, 'HGPS') > 0
%             myColorO = newcolors{5};
%             myTypeOuter = 'HGPS';
%         end
        currentLookedGroup = groupsCell(clusterLabels{myGroup}).RowNodeNames;
        aPol = [];
        aColors = [];
        aType = [];
        RiskCount = 0;
        MildCount = 0;
        SevereCount = 0;
        HCCount = 0;
        ProgeriaCount = 0;
        for lookedGroup=1:length(currentLookedGroup)
            currentWell = string(currentLookedGroup(lookedGroup));
%             if strfind(currentWell, 'Premanifest') > 0
%                 RiskCount = RiskCount + 1;
%                 myColor = newcolors{2};
%                 myType = 'Risk';
%             elseif strfind(currentWell, 'Mild') > 0
%                 MildCount = MildCount + 1;
%                 myColor = newcolors{3};
%                 myType = 'Mild';
%             elseif strfind(currentWell, 'Severe') > 0
%                 SevereCount = SevereCount + 1;
%                 myColor = newcolors{4};
%                 myType = 'Severe';
%             elseif strfind(currentWell, 'HC') > 0
%                 HCCount = HCCount + 1;
%                 myColor = newcolors{1};
%                 myType = 'HC';
%             elseif strfind(currentWell, 'HGPS') > 0
%                 ProgeriaCount = ProgeriaCount + 1;
%                 myColor = newcolors{5};
%                 myType = 'HGPS';
%             end
            number = str2double(regexp(currentWell,'\d*','Match', 'once'));
            aPol = [aPol one_dim_xys(number,:)'];
%             vColors = repmat({myColor},size(one_dim_xys_not_normalized(number),1),1);
%             aColors = [aColors; vColors];
%             vType = repmat({myType},size(one_dim_xys_not_normalized(number),1),1);
%             aType = [aType; vType];
        end
        if size(aPol,1) < numOfRandoms
            numOfRandomsIn = size(aPol,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
%         out = [];
%         indexBuild = 0;
%         typeOuter = 0;
%         typeElse = 0;
%         if strcmp(myTypeOuter, 'Risk')
%             usedOuter = RiskCount;
%             needed = RiskCount / max([MildCount SevereCount HCCount ProgeriaCount]);
%         elseif strcmp(myTypeOuter, 'Mild')
%             usedOuter = MildCount;
%             needed = MildCount / max([RiskCount SevereCount HCCount ProgeriaCount]);
%         elseif strcmp(myTypeOuter, 'Severe')
%             usedOuter = SevereCount;
%             needed = SevereCount / max([MildCount RiskCount HCCount ProgeriaCount]);
%         elseif strcmp(myTypeOuter, 'HC')
%             usedOuter = HCCount;
%             needed = HCCount / max([MildCount SevereCount RiskCount ProgeriaCount]);
%         elseif strcmp(myTypeOuter, 'HGPS')
%             usedOuter = ProgeriaCount;
%             needed = ProgeriaCount / max([MildCount SevereCount RiskCount HCCount]);
%         end
%         calc = 1 - 1/needed;
%         breakLoop = false;
%         while true
%             if breakLoop
%                 break;
%             end
%             if size(out, 2) == numOfRandomsIn
%                 breakLoop = true;
%                 break;
%             end
%             r = randi([1 size(aPol,1)],1,1);
%             if sum(ismember(out, r)) ~= 0
%                 continue;
%             end
%             if strcmp(aType(r),myTypeOuter)
%                
%                 if typeOuter > ceil(numOfRandomsIn*calc) + 2
%                     continue;
%                 end
%                 typeOuter = typeOuter + 1;
%                 out = [out r];
%             else
%                 if typeOuter > ceil(numOfRandomsIn*calc) + 2 || typeOuter == usedOuter
%                    typeElse = typeElse + 1;
%                     out = [out r]; 
%                 end
%             end
%             
%             
%         end
%         out = out(randperm(length(out)));
%         out = randperm(size(aPol,1),numOfRandomsIn);
%         outAll = randperm(size(aPol,1),size(aPol,1));
%         counter_x = 0;
%         counter_y = 0;
%         space = 500;
%         nexttile
%         hold on
        currNums = [];
        figure()
        mover = 0;
        for j=1:numOfRandomsIn
            currNum = round(size(aPol,1)*rand());
            while ismember(currNum, currNums) || currNum == 0
                currNum = ceil(size(aPol,1)*rand());
            end
            % getting sepcific cell trajectory
            xy=aPol(:,currNum); 
            xAxis = 0:10:10*73;
            c = linspace(1,10,length(xAxis));
            scatter(xAxis, xy + mover,sz, c, 'filled'); hold on;
            colorbar
            colormap jet
            plot(xAxis, xy + mover,'-','color','k', 'linewidth',0.5); hold on;
            title(['Activity Profile for ' currentGroupName])
            xlabel('Elapsed Time (750 minutes)') 
            ylabel('Z-Score Displacements') 
            set(gca, 'YTick', [])
%             xlim([-5 100])
%             col = aColors{currNum};
            mover = mover + 10;

            % zero cell position at first frame.
%             xy=xy-ones(size(xy(:,1)))*xy(1,:);
%             plot(xy(:,1)+counter_x*space,xy(:,2)+counter_y*space,'Color',col, 'LineWidth', 2);
%             counter_x = counter_x + 1;
%             if mod(counter_x,10)==0
%                 counter_x = 0;
%                 counter_y = counter_y + 1;
%             end
        end
%         set(gca,'xtick',[],'ytick',[]);
%         title(strcat("Trajectories for ", string(clusterNames{myGroup})));
%         axis square;
%         hold off
%         box on
%         ylim([-space space*5])
%         xlim([-space space*10])
%     % Show Sunplot of Cells
%         nexttile
%         cidk = jet(numOfRandomsIn);
%         for k=1:1:numOfRandomsIn
%              currNum = out(k);
%              xy=aPol{currNum};         
%              xy=xy-ones(size(xy(:,1)))*xy(1,:); % zero cell position at first frame.
%               xy1 = smooth(xy(:,1));
%              xy2 = smooth(xy(:,2));
%              plot(xy1,xy2,'-','color',[cidk(k,:) 1],'linewidth', 2); hold on; 
%         end
    end
    g = gcf;
    exportgraphics(g, strcat(folderCurrent, '\Ag_Activity_' + string(iBig) + '.png'), 'Resolution', 900)

    end
end
%% Lags and Trains
% binarize
% close all
meanArray = mean(one_dim_xys,2);
stdArray = std(one_dim_xys,0,2);
binarized = one_dim_xys>meanArray+0.25*stdArray;
TrainSum = [];
LagsSum = [];
TrainCV = [];
LagsCV = [];
grp1 = [];
labels = [];
if ~isempty(clusterLabels)
    
    for iBig=1:1:length(clusterLabels)
    starts = iBig;
    ends = iBig+0;
    if ends > length(clusterLabels)
        ends = length(clusterLabels); 
    end
    for myGroup=starts:ends
        if myGroup == 1
            myGroupUse = 4;
        elseif myGroup == 2
            myGroupUse = 2;
        elseif myGroup == 3
            myGroupUse = 1;
        elseif myGroup == 4
            myGroupUse = 3;
        end
        currentGroupName = clusterNames{myGroupUse};
        labels = [labels, currentGroupName];
        currentLookedGroup = groupsCell(clusterLabels{myGroupUse}).RowNodeNames;
        TrainsPol = [];
        LagsPol = [];
        for lookedGroup=1:length(currentLookedGroup)
            currentWell = string(currentLookedGroup(lookedGroup));
            number = str2double(regexp(currentWell,'\d*','Match', 'once'));
            TrainsPol = [TrainsPol; binarized(number,:)];
            LagsPol = [LagsPol; 1-binarized(number,:)];
        end
        TrainSum = [TrainSum, sum(TrainsPol, 2)'];
        LagsSum = [LagsSum, sum(LagsPol, 2)'];
        TrainCV = [TrainCV, (std(TrainsPol,0, 2)./mean(TrainsPol, 2))'];
        LagsCV = [LagsCV, (std(LagsPol,0, 2)./mean(LagsPol, 2))'];
%         if myGroup == 1
%             gr = (3)*ones(1,length(TrainsPol));
%         elseif myGroup == 2
%             gr = (0)*ones(1,length(TrainsPol));
%         elseif myGroup == 3
%             gr = (1)*ones(1,length(TrainsPol));
%         elseif myGroup == 4
%             gr = (2)*ones(1,length(TrainsPol));
%         end
        grp1 = [grp1, (myGroup-1)*ones(1,length(TrainsPol))];
%         TrainsAvg = mean(TrainsPol);
%         TrainsCV = std(TrainsPol); 
    end
    end
    figure()
    boxplot(TrainSum, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', labels)
    title('Boxplot for Trains');
    ylabel('Train Length');
    ylim([0 max(TrainSum)*1.005]);
    medians = grpstats(TrainSum, grp1, 'mean');
    x_locations = 1:numel(medians);
    hold on;
    plot(x_locations, medians, '-m');
    hold off;
    g = gcf;
    exportgraphics(g, strcat(folderCurrent, '\Ah_Trains_Length.png'), 'Resolution', 900)
    figure()
    boxplot(LagsSum, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', labels)
    title('Boxplot for Lags');
    ylabel('Lags Length');
    ylim([0 max(LagsSum)*1.05]);
    medians = grpstats(LagsSum, grp1, 'mean');
    x_locations = 1:numel(medians);
    hold on;
    plot(x_locations, medians, '-m');
    hold off;
    g = gcf;
    exportgraphics(g, strcat(folderCurrent, '\Ah_Lags_Length.png'), 'Resolution', 900)
    figure()
    boxplot(TrainCV, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', labels)
    title('Boxplot for CV Trains');
    ylabel('CV Train Length');
    ylim([0 max(TrainCV)*1.05]);
    medians = grpstats(TrainCV, grp1, 'mean');
    x_locations = 1:numel(medians);
    hold on;
    plot(x_locations, medians, '-m');
    hold off;
    g = gcf;
    exportgraphics(g, strcat(folderCurrent, '\Ah_CVTrains_Length.png'), 'Resolution', 900)
    figure()
    boxplot(LagsCV, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', labels)
    title('Boxplot for CV Lags');
    ylabel('CV Lags Length');
    ylim([0 max(LagsCV)*1.05]);
    medians = grpstats(LagsCV, grp1, 'mean');
    x_locations = 1:numel(medians);
    hold on;
    plot(x_locations, medians, '-m');
    hold off;
    g = gcf;
    exportgraphics(g, strcat(folderCurrent, '\Ah_CVLags_Length.png'), 'Resolution', 900)
end

%% Correlations between Motility and Activity
close all

figure() % by patient
uniqueCAPSHeatmap = uniqueCAPS;
uniqueCAPSHeatmap(1) = 'HC';
uniqueCAPSHeatmap(end) = 'HGPS';
YAxisNormed = bsxfun(@rdivide, YAxis, sum(YAxis, 1));
heatmap(uniqueCAPSHeatmap,XAxis, YAxisNormed, 'Colormap', parula);
%
title('Abundance of Activity cluster per Patient')
g = gcf;
saveas(g, strcat(folderCurrent, '\Aj_heatmap_patient.png'))
figure() % by motility clusters
M = zeros(length(GI_One_Dim), length(GI));
for iA=1:1:length(GI_One_Dim)
    currentClusterActivity = GI_One_Dim(iA).RowNodeNames;
    for iM=1:1:length(GI)
        currentClusterMotility = GI(iM).RowNodeNames;
        C = ismember(currentClusterActivity, currentClusterMotility);
        num_common = sum(C);
        M(iA, iM) = num_common/length(currentClusterActivity);% / length(currentClusterActivity);
    end
end
% heatmap(string(clusterNamesMotility),string(clusterNamesActivity), M/max(max(M)), 'Colormap', parula);
heatmap(string(clusterNamesMotility),string(clusterNamesActivity), M, 'Colormap', parula);
title('Abundance of Activity cluster per Motility cluster - Each row sums to 1.0')
g = gcf;
saveas(g, strcat(folderCurrent, '\Aj_heatmap_motility_Rows.png'))
figure()
M = zeros(length(GI_One_Dim), length(GI));
for iA=1:1:length(GI_One_Dim)
    currentClusterActivity = GI_One_Dim(iA).RowNodeNames;
    for iM=1:1:length(GI)
        currentClusterMotility = GI(iM).RowNodeNames;
        C = ismember(currentClusterMotility, currentClusterActivity);
        num_common = sum(C);
        M(iA, iM) = num_common/length(currentClusterMotility);% / length(currentClusterActivity);
    end
end
% heatmap(string(clusterNamesMotility),string(clusterNamesActivity), M/max(max(M)), 'Colormap', parula);
heatmap(string(clusterNamesMotility),string(clusterNamesActivity), M, 'Colormap', parula);
title('Abundance of Activity cluster per Motility cluster - Each column sums to 1.0')
g = gcf;
saveas(g, strcat(folderCurrent, '\Aj_heatmap_motility_Columns.png'))
%% by motility clusters per severence by CAP
close all
for iP=1:1:length(uniqueCAPS)
%     if uniqueCAPS(iP) ~= "126"
%        continue; 
%     end
    currentCAP = uniqueCAPS(iP) + "'";
    M = zeros(length(GI_One_Dim), length(GI));
    for iA=1:1:length(GI_One_Dim)
        currentClusterActivity = GI_One_Dim(iA).RowNodeNames;
        conntainsCAPActivity = contains(currentClusterActivity, currentCAP);
        currentClusterActivity = currentClusterActivity(conntainsCAPActivity);
        for iM=1:1:length(GI)
            currentClusterMotility = GI(iM).RowNodeNames;
            conntainsCAPMotility = contains(currentClusterMotility, currentCAP);
            currentClusterMotility = currentClusterMotility(conntainsCAPMotility);
            C = ismember(currentClusterActivity, currentClusterMotility);
            num_common = sum(C);
            if isempty(currentClusterMotility)
                if M(iA, iM) ~= 0
                    disp("ERROR");
                end
                M(iA, iM) = 0;
            else
                M(iA, iM) = num_common; % / length(currentClusterMotility);
            end
        end
    end
    figure()
    heatmap(string(clusterNamesMotility),string(clusterNamesActivity), M/max(max(M)), 'Colormap', parula);
    title(uniqueCAPS(iP));
    if iP == 1
        title('HC');
    end
    if iP == length(uniqueCAPS)
        title('HGPS');
    end
    g = gcf;
    saveas(g, strcat(folderCurrent, '\Ak_heatmap_CAP' + string(iP) + '.png'))
end
%% by motility clusters per severence
close all
for iP=1:1:length(uniqueCAPS)
    if str2double(uniqueCAPS(iP)) == 0 || str2double(uniqueCAPS(iP)) == 60 || str2double(uniqueCAPS(iP)) == 91 || str2double(uniqueCAPS(iP)) == 115 || str2double(uniqueCAPS(iP)) == 999
       M = zeros(length(GI_One_Dim), length(GI));
       counter = zeros(length(GI_One_Dim), length(GI));
    end
    currentCAP = uniqueCAPS(iP) + "'";

    for iM=1:1:length(GI)
        currentClusterMotility = GI(iM).RowNodeNames;
        conntainsCAPMotility = contains(currentClusterMotility, currentCAP);
        currentClusterMotility = currentClusterMotility(conntainsCAPMotility);
        for iA=1:1:length(GI_One_Dim)
            currentClusterActivity = GI_One_Dim(iA).RowNodeNames;
            conntainsCAPActivity = contains(currentClusterActivity, currentCAP);
            currentClusterActivity = currentClusterActivity(conntainsCAPActivity);
            C = ismember(currentClusterActivity, currentClusterMotility);
            num_common = sum(C);

            M(iA, iM) = M(iA, iM) + num_common;
            counter(iA, iM) = counter(iA, iM) + length(currentClusterMotility);
        end
            
    end
    if iP == length(uniqueCAPS)
        for iM=1:1:length(GI)
            for iA=1:1:length(GI_One_Dim)
                if counter(iA, iM) == 0
                    if M(iA, iM) ~= 0
                        disp("ERROR");
                    end
                    M(iA, iM) = 0;
                else
                    M(iA, iM) = M(iA, iM);% / counter(iA, iM);
                end
            end
        end
        figure()
        heatmap(string(clusterNamesActivity), string(clusterNamesMotility),M'/max(max(M')), 'Colormap', parula);
        title('HGPS');
        g = gcf;
        saveas(g, strcat(folderCurrent, '\Ak_heatmap_Severity' + string(iP) + '.png'))
    else
    if str2double(uniqueCAPS(iP+1)) == 0 || str2double(uniqueCAPS(iP+1)) == 60 || str2double(uniqueCAPS(iP+1)) == 91 || str2double(uniqueCAPS(iP+1)) == 115 || str2double(uniqueCAPS(iP+1)) == 999
        for iM=1:1:length(GI)
            for iA=1:1:length(GI_One_Dim)
                if counter(iA, iM) == 0
                    if M(iA, iM) ~= 0
                        disp("ERROR");
                    end
                    M(iA, iM) = 0;
                else
                    M(iA, iM) = M(iA, iM);% / counter(iA, iM);
                end
            end
        end
        figure()
        heatmap(string(clusterNamesActivity), string(clusterNamesMotility),M'/max(max(M')), 'Colormap', parula);
        title(uniqueCAPS(iP));
        if iP == 1
            title('HC');
        end
        if iP == length(uniqueCAPS)
            title('HGPS');
        end
        g = gcf;
        saveas(g, strcat(folderCurrent, '\Ak_heatmap_Severity' + string(iP) + '.png'))
    end
    end

end
%% Correlations between Activity and Patients
close all
% myCellTotal = [myCellHCCAPS; myCellRCAPS; myCellMCAPS; myCellSCAPS; myCellPCAPS];
% XAxis = categorical(myCellTotal(:,1));
% XAxis = reordercats(XAxis,flip(myCellTotal(:,1)));
% YAxis = str2double(myCellTotal(:,2:end));
