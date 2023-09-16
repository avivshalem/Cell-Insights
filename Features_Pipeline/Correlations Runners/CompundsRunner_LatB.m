%%% Correlations
%%%
load('\NatureAviv\Compunds2\LatB_LowHigh\xys_TRJs_MitoQ.mat')
load('\NatureAviv\Compunds2\LatB_LowHigh\xys_full_MitoQ.mat')
load('\NatureAviv\Compunds2\LatB_LowHigh\DataFileMitoQ.mat')

data = DataFileMitoQ;
data.CellType(data.CellType == "'LatBLow'") = "'LatB'";
data.CellType(data.CellType == "'LatBHigh'") = "'LatB'";

% find Wells
feature_vectors = [data.Dp, data.Dtot, data.Pp, data.Psi,...
    data.Pnp, data.Dnp, data.MSD10, ...
    data.Sp, data.Snp];

uniques = unique(data.CellType);
AllCountHCAll = find(data.CellType == uniques(1,:));
AllCountHCMediaAll = find(data.CellType == uniques(2,:));
AllCountHDAll = find(data.CellType == uniques(3,:));
AllCountHDMediaAll = find(data.CellType == uniques(4,:));
AllCountLatBAll = find(data.CellType == uniques(5,:));

%
% combinedIndex = combinedIndexCached;
% Correlations Loop Fixed X AXIS
folderCurrent = "Results\LatB";
mkdir(folderCurrent);

% feature_vectors = [data.Dp(combinedIndex), data.Dtot(combinedIndex), data.Pp(combinedIndex), data.Psi(combinedIndex),...
%     data.Pnp(combinedIndex), data.Dnp(combinedIndex), data.MSD10(combinedIndex), ...
%     data.Sp(combinedIndex), data.Snp(combinedIndex)];
% feature_names = [{'Dp'}, {'Dtot'},{'Pp'}, {'Psi'}, {'Pnp'}, {'Dnp'},...
%     {'MSD10'}, {'Sp'}, {'Snp'}];
% TF = isoutlier(feature_vectors, "quartiles", ThresholdFactor=1.5);
% TFRemoved = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
% % num = 8.8; 1
% % num = 6; 2
% % num = 5.5; 3
% % num = 2.75; 4 mean
num = 1.5;

FactorHC = num;
FactorPre = num;
FactorMild = num;
FactorSevere = num;
FactorHGPS = num;
% 
combinedIndex = ones(1,length(feature_vectors));
combinedIndexRemoved = true(size(combinedIndex));

% combinedIndexRemoved = false(size(combinedIndex));
% currentCombinedIndex = false(size(combinedIndex));
% currentCombinedIndex(AllCountHCAll) = combinedIndex(AllCountHCAll);
% TF = isoutlier(feature_vectors(currentCombinedIndex,:), "quartiles", ThresholdFactor=FactorHC);
% TFRemovedHC = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
% TFIndex = 1;
% for myIndex=1:length(currentCombinedIndex)
%     if currentCombinedIndex(myIndex)
%         if TFRemovedHC(TFIndex)
%             if combinedIndexRemoved(myIndex)
%                 disp("******** ERROR ***********");
%             end
%             combinedIndexRemoved(myIndex) = true;
%         end
%         TFIndex = TFIndex + 1;
%     end
% end
% currentCombinedIndex = false(size(combinedIndex));
% currentCombinedIndex(AllCountHDAll) = combinedIndex(AllCountHDAll);
% TF = isoutlier(feature_vectors(currentCombinedIndex,:), "quartiles", ThresholdFactor=FactorPre);
% TFRemovedPremanifest = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
% TFIndex = 1;
% for myIndex=1:length(currentCombinedIndex)
%     if currentCombinedIndex(myIndex)
%         if TFRemovedPremanifest(TFIndex)
%             if combinedIndexRemoved(myIndex)
%                 disp("******** ERROR ***********");
%             end
%             combinedIndexRemoved(myIndex) = true;
%         end
%         TFIndex = TFIndex + 1;
%     end
% end
% currentCombinedIndex = false(size(combinedIndex));
% currentCombinedIndex(AllCountMitoQAll) = combinedIndex(AllCountMitoQAll);
% TF = isoutlier(feature_vectors(currentCombinedIndex,:), "quartiles", ThresholdFactor=FactorMild);
% TFRemovedMild = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
% TFIndex = 1;
% for myIndex=1:length(currentCombinedIndex)
%     if currentCombinedIndex(myIndex)
%         if TFRemovedMild(TFIndex)
%             if combinedIndexRemoved(myIndex)
%                 disp("******** ERROR ***********");
%             end
%             combinedIndexRemoved(myIndex) = true;
%         end
%         TFIndex = TFIndex + 1;
%     end
% end
% currentCombinedIndex = false(size(combinedIndex));
% currentCombinedIndex(AllCountHCMediaAll) = combinedIndex(AllCountHCMediaAll);
% TF = isoutlier(feature_vectors(currentCombinedIndex,:), "quartiles", ThresholdFactor=FactorSevere);
% TFRemovedSevere = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
% TFIndex = 1;
% for myIndex=1:length(currentCombinedIndex)
%     if currentCombinedIndex(myIndex)
%         if TFRemovedSevere(TFIndex)
%             if combinedIndexRemoved(myIndex)
%                 disp("******** ERROR ***********");
%             end
%             combinedIndexRemoved(myIndex) = true;
%         end
%         TFIndex = TFIndex + 1;
%     end
% end
% currentCombinedIndex = false(size(combinedIndex));
% currentCombinedIndex(AllCountHDMediaAll) = combinedIndex(AllCountHDMediaAll);
% TF = isoutlier(feature_vectors(currentCombinedIndex,:), "quartiles", ThresholdFactor=FactorHGPS);
% TFRemovedHGPS = ~TF(:,1) & ~TF(:,2) & ~TF(:,3) & ~TF(:,4) & ~TF(:,5) & ~TF(:,6) & ~TF(:,7) & ~TF(:,8) & ~TF(:,9);
% TFIndex = 1;
% for myIndex=1:length(currentCombinedIndex)
%     if currentCombinedIndex(myIndex)
%         if TFRemovedHGPS(TFIndex)
%             if combinedIndexRemoved(myIndex)
%                 disp("******** ERROR ***********");
%             end
%             combinedIndexRemoved(myIndex) = true;
%         end
%         TFIndex = TFIndex + 1;
%     end
% end
combinedIndex = logical(combinedIndexRemoved);
combinedIndex(data.CellType == uniques(2,:)) = false;
combinedIndex(data.CellType == uniques(4,:)) = false;

% combinedIndexRemoved = [];
% TFIndex = 1;
% for myIndex=1:length(combinedIndex)
%     if combinedIndex(myIndex)
%         if TFRemoved(TFIndex)
%             combinedIndexRemoved = [combinedIndexRemoved; true];
%         else
%             combinedIndexRemoved = [combinedIndexRemoved; false];
%         end
%         TFIndex = TFIndex + 1;
%     else
%         combinedIndexRemoved = [combinedIndexRemoved; false];
%     end
% end
% combinedIndex = logical(combinedIndexRemoved);

for i=1:length(combinedIndex)
    if combinedIndex(i) == 0
        continue
    end
    if data.Psi(i) > 1000*median(data.Psi)
        combinedIndex(i) = 0;
    end
end
% AllCount1HC = find(data.CellType(combinedIndex) == uniques(1,:));
% AllCount1HGPS = find(data.CellType(combinedIndex) == uniques(2,:));
% AllCount1Mild = find(data.CellType(combinedIndex) == uniques(3,:));
% AllCount1Premanifest = find(data.CellType(combinedIndex) == uniques(4,:));
% AllCount1Severe = find(data.CellType(combinedIndex) == uniques(5,:));
%
feature_vectors = [data.Dp(combinedIndex), data.Dtot(combinedIndex), data.Pp(combinedIndex), data.Psi(combinedIndex),...
    data.Pnp(combinedIndex), data.Dnp(combinedIndex), data.MSD10(combinedIndex), ...
    data.Sp(combinedIndex), data.Snp(combinedIndex)];
feature_names = [{'Dp'}, {'Dtot'},{'Pp'}, {'Psi'}, {'Pnp'}, {'Dnp'},...
    {'MSD10'}, {'Sp'}, {'Snp'}];
% tempPatientWell = data.CAPPatient2(feature_vectors);
% new = [];
% uniquePatientsWell = unique(data.CAPPatient2(combinedIndex));
% for i=1:length(uniquePatientsWell)
%         currentData = find(tempPatientWell == uniquePatientsWell(i));
%         new = [new; mean(feature_vectors(currentData,:))];
% end
%
tempType = data.CellType(combinedIndex);
tempPatient = data.CAPPatient2(combinedIndex);
tempPatient(contains(string(tempPatient), "0'HCMedia")) = "0'HCMedia";
tempPatient(contains(string(tempPatient), "0'HC")) = "0'HC";

uniquePatients = unique(tempPatient);
uniquePatients = natsort(uniquePatients);

% DMSO
close all
tileFirst = 1 + floor(size(feature_vectors,2) / 2);

for iBig=1:4:size(feature_vectors,2)
starts = iBig;
ends = iBig+3;
if ends > size(feature_vectors,2)
    ends = size(feature_vectors,2); 
end
fig = figure('Position', get(0, 'Screensize'));
t = tiledlayout(2,2, "TileSpacing", "loose");
for f=starts:ends

    meanListHD = [];
    meanStdHD = [];
    meanListMQ = [];
    meanStdMQ = [];
    XTry = [];
    XMQTry = [];
    for i=1:length(uniquePatients)
        currentData = find(tempPatient == uniquePatients(i));
        currentFeatures = feature_vectors(currentData,f);
        currentType = unique(tempType(currentData));
        if length(currentType) > 1
            disp('ERROR');
        end
        if (contains(string(uniquePatients(i)), "'HD") || contains(string(uniquePatients(i)), "'HC"))...
                && ~contains(string(uniquePatients(i)), "Media")
            meanListHD(1+end) = mean(currentFeatures);
            meanStdHD(1+end) = std(currentFeatures);
            XTry(1+end) = str2double(regexp(string(uniquePatients(i)),'\d*','Match', 'once'));
        end
        if (contains(string(uniquePatients(i)), "'LatBHigh") || contains(string(uniquePatients(i)), "'HC"))...
                && ~contains(string(uniquePatients(i)), "Media")
            meanListMQ(1+end) = mean(currentFeatures);
            meanStdMQ(1+end) = std(currentFeatures);
            XMQTry(1+end) = str2double(regexp(string(uniquePatients(i)),'\d*','Match', 'once'));
        end
    end
    nexttile;
%     currX = 1:length(meanListHD);
%     XLabels = str2double(regexp(string(uniquePatients),'\d*','Match', 'once'));
%     XLabelsC = XLabels;
%     XLabels(XLabels==0) = 40;
%     XLabels(XLabels==999) = 161;
%     currX = [40 60 69 76 76 96 100 104 106 110 111 114 115 126 129 129 141];
%     currXMQ = [40 60 69 76 76 100 104 106 110 111 114 115 126 129 129 141];
    currX = sort(XTry);
    currX(1) = 40;
    currXMQ = sort(XMQTry);
    currXMQ(1) = 40;

    errorbar(currX, meanListHD,meanStdHD, 'bp', 'MarkerEdgeColor',[1 0 0],...
                'MarkerFaceColor',[0 .5 .5],...
                'LineWidth',0.25);
    hold on
    errorbar(currXMQ, meanListMQ,meanStdMQ, 'bp', 'MarkerEdgeColor',[0 1 0],...
                'MarkerFaceColor',[.5 0 .5],...
                'LineWidth',0.25);
    hold on
    CAPS = unique(currX);
    CAPS = CAPS(2:end)';
%     XT = [40, 60, 76, 94, 110, 119, 141, 161]; 
    XT = [40, CAPS'];
%     XLabels = ['HC', string(60), string(76), string(94), string(110), string(119), string(141), 'HGPS']; 
    XLabels = ['HC', string(CAPS')]; 
    set(gca, 'XTick', XT, 'XTickLabel', XLabels, 'fontweight','bold');
    pearson = corrcoef(currX,meanListHD);
    B = [currX(:) ones(size(currX(:)))] \ meanListHD(:);
    yfit = [currX(:) ones(size(currX(:)))]  * B;
    plot(currX, yfit, '-m')
    hold on
%     fitresult = fit(currX,meanListHD','exp1');
%     p11 = predint(fitresult,currX,0.8,'observation','off');
%     hold on;
%     plot(currX,p11,'k--')
%     hold off
%     xlim([min(XT)*0.8 max(XT)*1.05])
%     ylim([-0.1 1.15*max(meanStdHD + meanListHD)])
%     xlabel('CAP Score','fontweight','bold');
%     ylabel(string(feature_names(f)),'fontweight','bold');
    hold on
    currX = currXMQ;
    pearson2 = corrcoef(currX,meanListMQ);
    B = [currX(:) ones(size(currX(:)))] \ meanListMQ(:);
    yfit = [currX(:) ones(size(currX(:)))]  * B;
    plot(currX, yfit, '-c')
    hold on
%     fitresult = fit(currX,meanListMQ','exp1');
%     p11 = predint(fitresult,currX,0.8,'observation','off');
%     hold on;
%     plot(currX,p11,'k--')
    hold off
    xlim([min(XT)*0.8 max(XT)*1.05])
    ylim([-0.1 1.15*max(meanStdHD + meanListHD)])
    xlabel('CAP Score','fontweight','bold');
    ylabel(string(feature_names(f)),'fontweight','bold');
    if length(pearson) == 1
        title(strcat('HD Pearson Coefficient is', {' '}, num2str(pearson), {' '}, 'LatB Pearson Coefficient is', {' '}, num2str(pearson2)),'fontweight','bold');
    else
        title(strcat('HD Pearson Coefficient is', {' '}, num2str(pearson(2)), {' '}, 'LatB Pearson Coefficient is', {' '}, num2str(pearson2(2))),'fontweight','bold');
    end

end
F    = getframe(fig);
exportgraphics(t,strcat(folderCurrent, '\_z_CAP_Correlations_', num2str(iBig), '.png'),'Resolution',500)
end
% %% Media
% close all
% tileFirst = 1 + floor(size(feature_vectors,2) / 2);
% 
% for iBig=1:4:size(feature_vectors,2)
% starts = iBig;
% ends = iBig+3;
% if ends > size(feature_vectors,2)
%     ends = size(feature_vectors,2); 
% end
% fig = figure('Position', get(0, 'Screensize'));
% t = tiledlayout(2,2, "TileSpacing", "loose");
% for f=starts:ends
% 
%     meanListHD = [];
%     meanStdHD = [];
%     meanListMQ = [];
%     meanStdMQ = [];
%     XTry = [];
%     XMQTry = [];
%     for i=1:length(uniquePatients)
%         currentData = find(tempPatient == uniquePatients(i));
%         currentFeatures = feature_vectors(currentData,f);
%         currentType = unique(tempType(currentData));
%         if length(currentType) > 1
%             disp('ERROR');
%         end
%         if (contains(string(uniquePatients(i)), "HDMedia") || contains(string(uniquePatients(i)), "HCMedia"))...
%                 && contains(string(uniquePatients(i)), "Media")
%             meanListHD(1+end) = mean(currentFeatures);
%             meanStdHD(1+end) = std(currentFeatures);
%             XTry(1+end) = str2double(regexp(string(uniquePatients(i)),'\d*','Match', 'once'));
% 
%         end
%         if (contains(string(uniquePatients(i)), "'MitoQ") || contains(string(uniquePatients(i)), "HCMedia"))...
%                 && contains(string(uniquePatients(i)), "Media")
%             meanListMQ(1+end) = mean(currentFeatures);
%             meanStdMQ(1+end) = std(currentFeatures);
%             XMQTry(1+end) = str2double(regexp(string(uniquePatients(i)),'\d*','Match', 'once'));
% 
%         end
%     end
%     nexttile;
%     currX = 1:length(meanListHD);
%     XLabels = str2double(regexp(string(uniquePatients),'\d*','Match', 'once'));
%     XLabelsC = XLabels;
%     XLabels(XLabels==0) = 40;
% %     XLabels(XLabels==999) = 161;
%     currX = sort(XTry);
%     currX(1) = 40;
%     currXMQ = sort(XMQTry);
%     currXMQ(1) = 40;
% 
%     errorbar(currX, meanListHD,meanStdHD, 'bp', 'MarkerEdgeColor',[1 0 0],...
%                 'MarkerFaceColor',[0 .5 .5],...
%                 'LineWidth',0.25);
%     hold on
%     errorbar(currXMQ, meanListMQ,meanStdMQ, 'bp', 'MarkerEdgeColor',[0 1 0],...
%                 'MarkerFaceColor',[.5 0 .5],...
%                 'LineWidth',0.25);
%     hold on
%     CAPS = unique(currX);
%     CAPS = CAPS(2:end)';
% %     XT = [40, 60, 76, 94, 110, 119, 141, 161]; 
%     XT = [40, CAPS];
% %     XLabels = ['HC', string(60), string(76), string(94), string(110), string(119), string(141), 'HGPS']; 
%     XLabels = ['HC', string(CAPS)]; 
%     set(gca, 'XTick', XT, 'XTickLabel', XLabels, 'fontweight','bold');
%     pearson = corrcoef(currX,meanListHD);
%     B = [currX(:) ones(size(currX(:)))] \ meanListHD(:);
%     yfit = [currX(:) ones(size(currX(:)))]  * B;
%     plot(currX, yfit, '-m')
%     hold on
% %     fitresult = fit(currX,meanListHD','exp1');
% %     p11 = predint(fitresult,currX,0.8,'observation','off');
% %     hold on;
% %     plot(currX,p11,'k--')
% %     hold off
% %     xlim([min(XT)*0.8 max(XT)*1.05])
% %     ylim([-0.1 1.15*max(meanStdHD + meanListHD)])
% %     xlabel('CAP Score','fontweight','bold');
% %     ylabel(string(feature_names(f)),'fontweight','bold');
%     hold on
%     currX = currXMQ;
%     pearson2 = corrcoef(currX,meanListMQ);
%     B = [currX(:) ones(size(currX(:)))] \ meanListMQ(:);
%     yfit = [currX(:) ones(size(currX(:)))]  * B;
%     plot(currX, yfit, '-c')
%     hold on
% %     fitresult = fit(currX,meanListMQ','exp1');
% %     p11 = predint(fitresult,currX,0.8,'observation','off');
% %     hold on;
% %     plot(currX,p11,'k--')
%     hold off
%     xlim([min(XT)*0.8 max(XT)*1.05])
%     ylim([-0.1 1.15*max(meanStdHD + meanListHD)])
%     xlabel('CAP Score','fontweight','bold');
%     ylabel(string(feature_names(f)),'fontweight','bold');
%     if length(pearson) == 1
%         title(strcat('HD Pearson Coefficient is', {' '}, num2str(pearson), {' '}, 'MitoQ Pearson Coefficient is', {' '}, num2str(pearson2)),'fontweight','bold');
%     else
%         title(strcat('HD Pearson Coefficient is', {' '}, num2str(pearson(2)), {' '}, 'MitoQ Pearson Coefficient is', {' '}, num2str(pearson2(2))),'fontweight','bold');
%     end
% 
% end
% F    = getframe(fig);
% exportgraphics(t,strcat(folderCurrent, '\_z_CAP_Correlations_Media', num2str(iBig), '.png'),'Resolution',500)
% end
% %% Both
% close all
% tileFirst = 1 + floor(size(feature_vectors,2) / 2);
% 
% for iBig=1:4:size(feature_vectors,2)
% starts = iBig;
% ends = iBig+3;
% if ends > size(feature_vectors,2)
%     ends = size(feature_vectors,2); 
% end
% fig = figure('Position', get(0, 'Screensize'));
% t = tiledlayout(2,2, "TileSpacing", "loose");
% for f=starts:ends
% 
%     meanListHD = [];
%     meanStdHD = [];
%     meanListMQ = [];
%     meanStdMQ = [];
%     meanListHDMedia = [];
%     meanStdHDMedia = [];
%     meanListMQMedia = [];
%     meanStdMQMedia = [];
%     for i=1:length(uniquePatients)
%         currentData = find(tempPatient == uniquePatients(i));
%         currentFeatures = feature_vectors(currentData,f);
%         currentType = unique(tempType(currentData));
%         if length(currentType) > 1
%             disp('ERROR');
%         end
%         if contains(string(uniquePatients(i)), "HDMedia") || contains(string(uniquePatients(i)), "HCMedia")
%             meanListHDMedia(1+end) = mean(currentFeatures);
%             meanStdHDMedia(1+end) = std(currentFeatures);
%         end
%         if contains(string(uniquePatients(i)), "'MitoQ") || contains(string(uniquePatients(i)), "HCMedia")
%             meanListMQMedia(1+end) = mean(currentFeatures);
%             meanStdMQMedia(1+end) = std(currentFeatures);
%         end
%         if contains(string(uniquePatients(i)), "'HD") || contains(string(uniquePatients(i)), "'HC")
%             meanListHD(1+end) = mean(currentFeatures);
%             meanStdHD(1+end) = std(currentFeatures);
%         end
%         if contains(string(uniquePatients(i)), "'MitoQ") || contains(string(uniquePatients(i)), "'HC")
%             meanListMQ(1+end) = mean(currentFeatures);
%             meanStdMQ(1+end) = std(currentFeatures);
%         end
%     end
%     nexttile;
%     currX = 1:length(meanListHDMedia);
%     XLabels = str2double(regexp(string(uniquePatients),'\d*','Match', 'once'));
%     XLabelsC = XLabels;
%     XLabels(XLabels==0) = 40;
% %     XLabels(XLabels==999) = 161;
%     currXMedia = [40 96 100 110 111];
%     currXMQ = [40 60 69 76 76 100 104 106 110 111 114 115 126 129 129 141];
% 
%     errorbar(currXMedia, meanListHDMedia,meanStdHDMedia, 'bp', 'MarkerEdgeColor',[1 0 0],...
%                 'MarkerFaceColor',[0 .5 .5],...
%                 'LineWidth',0.25);
%     hold on
%     errorbar(currXMQ, meanListMQMedia,meanStdMQMedia, 'bp', 'MarkerEdgeColor',[0 1 0],...
%                 'MarkerFaceColor',[.5 0 .5],...
%                 'LineWidth',0.25);
%     hold on
%     currXDMSO = [40 60 69 76 76 96 100 104 106 110 111 114 115 126 129 129 141];
%     errorbar(currXDMSO, meanListHD,meanStdHD, 'bp', 'MarkerEdgeColor',[0 0 1],...
%                 'MarkerFaceColor',[.5 .5 0],...
%                 'LineWidth',0.25);
%     hold on
%     errorbar(currXMQ, meanListMQ,meanStdMQ, 'bp', 'MarkerEdgeColor',[1 1 1],...
%                 'MarkerFaceColor',[.5 .5 .5],...
%                 'LineWidth',0.25);
%     hold on
%     CAPS = unique(XLabelsC);
%     CAPS = CAPS(2:end)';
% %     XT = [40, 60, 76, 94, 110, 119, 141, 161]; 
%     XT = [40, CAPS];
% %     XLabels = ['HC', string(60), string(76), string(94), string(110), string(119), string(141), 'HGPS']; 
%     XLabels = ['HC', string(CAPS)]; 
%     currX = currXMedia;
%     set(gca, 'XTick', XT, 'XTickLabel', XLabels, 'fontweight','bold');
%     pearson = corrcoef(currX,meanListHDMedia);
%     B = [currX(:) ones(size(currX(:)))] \ meanListHDMedia(:);
%     yfit = [currX(:) ones(size(currX(:)))]  * B;
%     plot(currX, yfit, '-m')
%     hold on
% %     fitresult = fit(currX,meanListHD','exp1');
% %     p11 = predint(fitresult,currX,0.8,'observation','off');
% %     hold on;
% %     plot(currX,p11,'k--')
% %     hold off
% %     xlim([min(XT)*0.8 max(XT)*1.05])
% %     ylim([-0.1 1.15*max(meanStdHD + meanListHD)])
% %     xlabel('CAP Score','fontweight','bold');
% %     ylabel(string(feature_names(f)),'fontweight','bold');
%     hold on
%     currX = currXMQ;
%     pearson2 = corrcoef(currX,meanListMQMedia);
%     B = [currX(:) ones(size(currX(:)))] \ meanListMQMedia(:);
%     yfit = [currX(:) ones(size(currX(:)))]  * B;
%     plot(currX, yfit, '-c')
%     hold on
% %     fitresult = fit(currX,meanListMQ','exp1');
% %     p11 = predint(fitresult,currX,0.8,'observation','off');
% %     hold on;
% %     plot(currX,p11,'k--')
%     currX = currXDMSO;
%     pearson3 = corrcoef(currX,meanListHD);
%     B = [currX(:) ones(size(currX(:)))] \ meanListHD(:);
%     yfit = [currX(:) ones(size(currX(:)))]  * B;
%     plot(currX, yfit, '-r')
%     hold on
% %     fitresult = fit(currX,meanListHD','exp1');
% %     p11 = predint(fitresult,currX,0.8,'observation','off');
% %     hold on;
% %     plot(currX,p11,'k--')
% %     hold off
% %     xlim([min(XT)*0.8 max(XT)*1.05])
% %     ylim([-0.1 1.15*max(meanStdHD + meanListHD)])
% %     xlabel('CAP Score','fontweight','bold');
% %     ylabel(string(feature_names(f)),'fontweight','bold');
%     hold on
%     currX = currXMQ;
%     pearson4 = corrcoef(currX,meanListMQ);
%     B = [currX(:) ones(size(currX(:)))] \ meanListMQ(:);
%     yfit = [currX(:) ones(size(currX(:)))]  * B;
%     plot(currX, yfit, '-y')
%     hold on
% %     fitresult = fit(currX,meanListMQ','exp1');
% %     p11 = predint(fitresult,currX,0.8,'observation','off');
% %     hold on;
% %     plot(currX,p11,'k--')
%     hold off
%     xlim([min(XT)*0.8 max(XT)*1.05])
%     ylim([-0.1 1.15*max(meanStdHD + meanListHD)])
%     xlabel('CAP Score','fontweight','bold');
%     ylabel(string(feature_names(f)),'fontweight','bold');
%     if length(pearson) == 1
%         title(strcat('HD Media Pearson Coefficient is', {' '}, num2str(pearson), {' '}, 'MitoQ Media Pearson Coefficient is', {' '}, num2str(pearson2), {' '}, ...
%             'HD DMSO Pearson Coefficient is', {' '}, num2str(pearson3), {' '}, 'MitoQ DMSO Pearson Coefficient is', {' '}, num2str(pearson4)),'fontweight','bold');
%     else
%         title(strcat('HD Pearson Coefficient is', {' '}, num2str(pearson(2)), {' '}, 'MitoQ Pearson Coefficient is', {' '}, num2str(pearson2(2)), {' '}, ...
%             'HD Pearson Coefficient is', {' '}, num2str(pearson3(2)), {' '}, 'MitoQ Pearson Coefficient is', {' '}, num2str(pearson4(2))),'fontweight','bold');
%     end
% 
% end
% F    = getframe(fig);
% exportgraphics(t,strcat(folderCurrent, '\_z_CAP_Correlations_', num2str(iBig), '.png'),'Resolution',500)
% end
%%
% % close all
% %
% % xysMain=get_trajfile; %(DATA.mat)
% % xysMain_progeria=get_trajfile; %(DATA.xlsx)
% 
% %%%
% idxComb = combinedIndex;
xysMain = xys_TRJs_MitoQ(combinedIndex);
xys_full = xys_full_MitoQ(combinedIndex);
% xys = xysMain(idxComb);

% newcolors = {'#4285F4','#DB4437','#F4B400','#DB4437' '#330072', '#016773'};
newcolors = {'#4285F4','#0F9D58','#F4B400','#DB4437' '#330072', '#016773'};

newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
% newcolorsRGB = {[66 133 244]/255,[219 68 55]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};


xys = cell(length(xys_full),1);
for i=1:length(xys_full)
    current = xys_full{i};
    xys{i} = xys_full{i}(:,[3,4]);
end

% MSD vs TLAG
% close all force
param.showfig=1;
param.saveres=0;
param.markertype='.-';
param.outfigurenum=301;
param.dim=2;
param.linear=0;
param.aviv=0;

temp_cellType = data.CellType(combinedIndex);

a1Pol = xys(temp_cellType == "'HC'");
a2Pol = xys(temp_cellType == "'HD'");
a3Pol = xys(temp_cellType == "'LatB'");
% a5Pol = xys(temp_cellType == "HDMedia'");

get_MSD(a1Pol, 10, param);
hold on;
param.markertype='.-';
get_MSD(a2Pol, 10, param);
hold on;
param.markertype='.-'; 
get_MSD(a3Pol, 10, param);
hold on;
% param.markertype='.-'; 
% get_MSD(a4Pol, 10, param);
% hold on;
% param.markertype='.-';
% get_MSD(a4Pol, 10, param);
% hold on;
% param.markertype='.-'; 
% get_MSD(a5Pol, 10, param);
% hold on;

% colororder(newcolors)
param.markertype='k:';
param.linear=1;
get_MSD(a1Pol, 10, param);
hold off;
title('Mean Squared Displacement');

ylabel('MSD (µm2)') ;
xlabel('Time lag (min)');

% legend('HC','HD' ,'MitoQ', 'HCMedia', 'HDMedia', 'Location','bestoutside');
legend('HC','HD' ,'LatB 60nM','Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';

saveas(g, strcat(folderCurrent, '\z_MSD.png'))
%%
%%% MSD vs TLAG PATIENT
% close all force
figure
param.showfig=1;
param.saveres=0;
param.aviv=1;
param.outfigurenum=301;
param.dim=2;
param.linear=0;

temp_Patient = data.CAP;
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
legend(name, 'Location', 'bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_MSDPatients.png'))
%%
% dTheta PDF
close all
a1Pol = xys(temp_cellType == "'HC'");
n = size(a1Pol, 1);  % number of elements in a1Pol
perm_indices = randperm(n);
half_indices = perm_indices(1:round(n / 10));
a1Pol = a1Pol(half_indices);

a2Pol = xys(temp_cellType == "'HD'");
n = size(a2Pol, 1);  % number of elements in a1Pol
perm_indices = randperm(n);
half_indices = perm_indices(1:round(n / 10));
a2Pol = a2Pol(half_indices);

a3Pol = xys(temp_cellType == "'LatB'");
n = size(a3Pol, 1);  % number of elements in a1Pol
perm_indices = randperm(n);
half_indices = perm_indices(1:round(n / 2));
a3Pol = a3Pol(half_indices);

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
saveas(g, strcat(folderCurrent, '\mqz_Theta_PDF_HC.png'))
% hold on;
% param.markertype=newcolorsRGB{2};
% figure();
param.outfigurenum = param.outfigurenum + 1;
get_dtheta_PDF(a2Pol, tloi, param);
% legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
title('HD dθ PDF');
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
saveas(g, strcat(folderCurrent, '\mqz_Theta_PDF_HD.png'))
% hold on;
% param.markertype=newcolorsRGB{3};
% figure();
param.outfigurenum = param.outfigurenum + 1;
get_dtheta_PDF(a3Pol, tloi, param);
% legend({'10 Minutes', '50 Minutes', '100 Minutes', '150 Minutes', '200 Minutes', '250 Minutes', '300 Minutes', '350 Minutes', '400 Minutes', '450 Minutes', });
title('LatB  60nM dθ PDF');
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
saveas(g, strcat(folderCurrent, '\mqz_Theta_PDF_LatB.png'))

%%
%%% dR Polarity
close all
param.showfig=1;
param.saveres=0;
param.dim=2;
param.outfigurenum=305;
param.binnum=35;

param.markertype=newcolorsRGB{1};
get_dR_polarity(a1Pol, 1, param);
hold on;

param.markertype=newcolorsRGB{3};
get_dR_polarity(a2Pol, 1, param);
hold on;
param.markertype=newcolorsRGB{4};
get_dR_polarity(a3Pol, 1, param);
hold on;
% param.markertype=newcolorsRGB{2};
% get_dR_polarity(a4Pol, 1, param);
% hold on;
% param.markertype=newcolorsRGB{5};
% get_dR_polarity(a5Pol, 1, param);
% hold on;
% colororder(newcolors)

title('dR Polarity');
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_Theta.png'))
%%
% legend('HC','Premanifest' ,'Mild', 'Severe','Progeria','Location','bestoutside');
%
%%% dR Polarity Patient
close all
param.showfig=1;
param.saveres=0;
param.dim=2;
param.outfigurenum=305;
param.binnum=35;


temp_Patient = data.CAP;
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
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_ThetaPatients.png'))
% legend(name, 'Location', 'bestoutside');
% legend('HC','Premanifest' ,'Mild', 'Severe','Location','bestoutside');
%%
%%% ACF vs TLAG
close all force
param.dim=2; 
param.tlag=1;
param.saveres=0;
param.showfig=1;       
param.markertype=''; % HC
param.outfigurenum=303;   
param.aviv=1;
param.aviv3 = 0;

get_ACF1(a1Pol, 10, param);
hold on;
param.markertype=''; % HD
get_ACF1(a2Pol, 10, param);
hold on;
param.markertype=''; % HD
get_ACF1(a3Pol, 10, param);
hold on;
% param.markertype=''; % HD
% get_ACF1(a4Pol, 10, param);
% hold on;
% param.markertype=''; % HD
% get_ACF1(a5Pol, 10, param);
% hold on;

% colororder(newcolors)

title('Auto-Correlation Function');
ylabel('ACF') ;
xlabel('Time lag (min)');
legend('HC','HD' ,'LatB  60nM','Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_ACF.png'))
%%
%%% ACF vs TLAG PATIENT
close all force
param.dim=2; 
param.tlag=1;
param.saveres=0;
param.showfig=1;       
param.markertype='m-'; % HC
param.outfigurenum=303;   
param.aviv=0;


temp_Patient = data.CAP;
cidk = jet(length(unique(temp_Patient)));
name = sort(unique(temp_Patient), 'descend');
for i=1:length(unique(temp_Patient))
    
    temp_pat = i;

    clust_idx = find(temp_Patient == name(i));
    check = find(temp_Patient == name(temp_pat));
 
    aPol = xys(check);
    
    param.aviv2 = 0;
    param.markertype=cidk(length(unique(temp_Patient)) + 1 - i,:);
    get_ACF1(aPol, 10, param);
    hold on;

    
end

title('Auto-Correlation Function by CAP Scores');
ylabel('ACF') ;
xlabel('Time lag (min)');

name = string(name);
legend(name, 'Location', 'bestoutside');
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_ACFPatients.png'))
%%
%%% dR PDF
close all force
param.dxmax=50;
param.binn=70;            
param.showfig=1;
param.saveres=0;
param.dim=2;            
param.outfigurenum=302;
param.markertype='o'; % HC
param.aviv=0;
param.aviv3=1;

param.markertype=':';
get_dR_PDF(a1Pol, 10, param);
hold on;
param.markertype=':'; % HD
get_dR_PDF(a2Pol, 10, param);
hold on;
param.markertype=':'; % HD
get_dR_PDF(a3Pol, 10, param);
hold on;
% param.markertype=':'; % HD
% get_dR_PDF(a4Pol, 10, param);
% hold on;
% param.markertype=':'; % HD
% get_dR_PDF(a5Pol, 10, param);
% hold on;
% colororder(newcolors)

title('PDF Cellular Displacements');
ylabel('Occurrence') ;
xlabel('Displacement (µm)');
legend('HC','HD' ,'LatB  60nM','Location','bestoutside');
    axis square;
    box on
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_PDF.png'))
%%
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

temp_Patient = data.CAP;
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
name(1) = 'HD - 110.9399';
name(2) = 'HD - 110.0154';
name(3) = 'HD - 100.1541';
name(4) = 'HD - 96.1479';
name(5) = 'MitoQ - 110.9399';
name(6) = 'MitoQ - 110.0154';
name(7) = 'MitoQ - 100.1541';
name(end) = 'HC';
legend(name, 'Location', 'bestoutside');
g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\z_PDFPatients.png'))
%%
xysMain = xys_TRJs_MitoQ(combinedIndex);
xys_full = xys_full_MitoQ(combinedIndex);
xys = cell(length(xys_full),1);
for i=1:length(xys_full)
    current = xys_full{i};
    xys{i} = xys_full{i}(:,[3,4]);
end
% xys = xys;
xys = xysMain;

temp_cellType = data.CellType(combinedIndex);

a1Pol = xys(temp_cellType == "'HC'");
a2Pol = xys(temp_cellType == "'HD'");
a3Pol = xys(temp_cellType == "'LatB'");

% let's assume classes are stored in a cell array
classes = {a1Pol, a2Pol, a3Pol};
classNames = {'HC', 'HD', 'LatB'};

% initialize empty arrays to store all radius of gyrations and corresponding labels
allRadiusOfGyration = [];
allLabels = [];

for i = 1:length(classes)
    % get the trajectories for this class
    trajectories = classes{i};
    
    % initialize an array to store the radius of gyration for this class
    radiusOfGyration = zeros(length(trajectories), 1);
    
    for j = 1:length(trajectories)
        trajectory = trajectories{j};
        x = trajectory(:, 1) - mean(trajectory(:, 1));
        y = trajectory(:, 2) - mean(trajectory(:, 2));
        
        % compute the radius of gyration
        radiusOfGyration(j) = sqrt(mean(x.^2 + y.^2));
    end
    
    % add the radius of gyration and labels for this class to the overall arrays
    allRadiusOfGyration = [allRadiusOfGyration; radiusOfGyration];
    allLabels = [allLabels; repmat(classNames(i), length(radiusOfGyration), 1)];
end

% now, you can create a boxplot for each class
figure
boxplot(allRadiusOfGyration, allLabels);
title('Radius of gyration for each class');

%%
%%% Draw TRJ

normalized_features = log(feature_vectors);
tempPatientWell = data.CAPPatient2(combinedIndex);
new = [];
uniquePatientsWell = unique(data.CAPPatient2(combinedIndex));
for i=1:length(uniquePatientsWell)
        currentData = find(tempPatientWell == uniquePatientsWell(i));
        new = [new; mean(normalized_features(currentData,:))];
end
z = zscore(normalized_features);
zWell = zscore(new);
%
xys = xysMain;
a1Pol = xys(temp_cellType == "'HC'");
a2Pol = xys(temp_cellType == "'HD'");
a3Pol = xys(temp_cellType == "'LatB'");
% a4Pol = xys(temp_cellType == "'Premanifest'");
% a5Pol = xys(temp_cellType == "'HGPS'");
%%% Show Random Cells
close all force
numOfRandoms = 50;

for SMR=1:5
    if SMR == 5 || SMR == 3
        continue;
    end
    fig = figure('Position', get(0, 'Screensize'));
    t = tiledlayout(3,1, 'TileSpacing','tight');
    if SMR == 1
        cells = a2Pol;
        titleSMR = {"Trajectories", "HD"};
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
        titleSMR = {"Trajectories", "LatB  60nM"};
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
        titleSMR = {"Trajectories", "LatBHigh"};
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
%     elseif SMR == 5
%         cells = a5Pol;
%         titleSMR = strcat("HGPS Trajectories");
%         color = newcolors{5};
%         if size(cells,1) < numOfRandoms
%             numOfRandomsIn = size(cells,1);
%         else
%             numOfRandomsIn = numOfRandoms;
%         end
%         out5 = randperm(size(cells,1),numOfRandomsIn);
%         out5All = randperm(size(cells,1),size(cells,1));
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
%         elseif SMR == 5
%             currNum = out5(j);
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
        titleSMR = strcat("HD Sunplot Trajectories");
        color = newcolors{4};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 2
        cells = a3Pol;
        titleSMR = strcat("LatBLow Sunplot Trajectories");
        color = newcolors{3};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 3
        cells = a4Pol;
        titleSMR = strcat("LatBHigh Sunplot Trajectories");
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
%     elseif SMR == 5
%         cells = a5Pol;
%         titleSMR = strcat("HGPS Sunplot Trajectories");
%         color = newcolors{5};
%         if size(cells,1) < numOfRandoms
%             numOfRandomsIn = size(cells,1);
%         else
%             numOfRandomsIn = numOfRandoms;
%         end
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
%         elseif SMR == 5
%             currNum = out5(k);
        end
         xy=cells{currNum};         
         xy=xy-ones(size(xy(:,1)))*xy(1,:); % zero cell position at first frame.
         xy1 = smooth(xy(:,1));
         xy2 = smooth(xy(:,2));
         plot(xy1,xy2,'-','color',[cidk(k,:) 1],'linewidth',2); hold on; 
    end
%     title(titleSMR);
    axis square
        box on
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-350 350])
    xlim([-350 350])
    text(90, 300, strcat("N = ", string(length(cells))),'FontSize', 16, 'FontWeight', 'bold');

    % New Plt medged

    if SMR == 1
        cells = a2Pol;
        titleSMR = strcat("HD Merged Trajectory");
        color = newcolors{4};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 2
        cells = a3Pol;
        titleSMR = strcat("LatBLow Merged Trajectory");
        color = newcolors{3};
        if size(cells,1) < numOfRandoms
            numOfRandomsIn = size(cells,1);
        else
            numOfRandomsIn = numOfRandoms;
        end
    elseif SMR == 3
        cells = a4Pol;
        titleSMR = strcat("LatBHigh Merged Trajectory");
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
%     elseif SMR == 5
%         cells = a5Pol;
%         titleSMR = strcat("HGPS Merged Trajectory");
%         color = newcolors{5};
%         if size(cells,1) < numOfRandoms
%             numOfRandomsIn = size(cells,1);
%         else
%             numOfRandomsIn = numOfRandoms;
%         end
    end
    nexttile
    cidk = jet(numOfRandomsIn);
    xy = 0;
    for k=1:1:numOfRandomsIn
        if SMR == 1
            currNum = out1All(k);
        elseif SMR == 2
            currNum = out2All(k);
        elseif SMR == 3
            currNum = out3All(k);
        elseif SMR == 4
            currNum = out4All(k);
%         elseif SMR == 5
%             currNum = out5All(k);
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
        box on
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-3 3])
    xlim([-3 3])
    exportgraphics(t,strcat(folderCurrent, 'a_', string(SMR),'_TRJFirst.jpg'),'Resolution',500)

end

%% Activity for Groups
%
% Prepare 1-dimensional displacement trajectories

xys = xysMain;

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

close all
numOfRandoms = 10;
sz = 10;

a1Pol_Activity = one_dim_xys(temp_cellType == "'HC'",:);
a2Pol_Activity = one_dim_xys(temp_cellType == "'HD'",:);
a3Pol_Activity = one_dim_xys(temp_cellType == "'LatB'",:);

classes = {a1Pol_Activity, a2Pol_Activity, a3Pol_Activity};
classNames = {'HC', 'HD', 'LatB  60nM'};
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
a2Pol_Activity = binarized(temp_cellType == "'HD'",:);
a3Pol_Activity = binarized(temp_cellType == "'LatB'",:);

classes = {a1Pol_Activity, a2Pol_Activity, a3Pol_Activity};
classNames = {'HC', 'HD', 'LatB  60nM'};

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
boxplot(TrainSum, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', classNames)
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

%%
close all force
%%% Show Random Cells Per CAP
numOfRandoms = 50;

temp_Patient = data.CAP;
name = sort(unique(temp_Patient), 'descend');
for i=1:length(unique(temp_Patient))
    fig = figure('Position', get(0, 'Screensize'));
    t = tiledlayout(3,1, 'TileSpacing','Compact');
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
        title(strcat("Trajectories for HC Cells"));
    elseif stringName == 999
        title(strcat("Trajectories for HGPS Cells"));
    else
        title(strcat("Trajectories for CAP Score: ", string(name(i))));
    end
    axis square;
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
    if stringName == 0
        title(strcat("Sunplot Trajectories for HC Cells"));
    elseif stringName == 999
        title(strcat("Sunplot Trajectories for HGPS Cells"));
    else
        title(strcat("Sunplot Trajectories for CAP Score: ", string(name(i))));
    end
    axis square
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-350 350])
    xlim([-350 350])
    
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
    if stringName == 0
        title(strcat("Merged Trajectory for HC Cells"));
    elseif stringName == 999
        title(strcat("Merged Trajectory for HGPS Cells"));
    else
        title(strcat("Merged Trajectory for CAP Score: ", string(name(i))));
    end
    axis square
    hold off;
    set(gca,'xtick',[],'ytick',[]);
    ylim([-30 30])
    xlim([-30 30])
    exportgraphics(t,strcat(folderCurrent, 'b_', string(i),'_TRJCAP.jpg'),'Resolution',500)

end
%
close all force
%% 2d look at features by Cells
close all
clc
num_dims = 2; % Reduce the data to 2 dimensions
% labels = regexp(string(uniquePatientsWell),'\d*','Match', 'once');
tempPatient(contains(string(tempPatient),"'HDMedia")) = 'HDMedia';
tempPatient(contains(string(tempPatient),"0'HCMedia")) = 'HCMedia';
tempPatient(contains(string(tempPatient),'LatBLow')) = 'LatBLow';
tempPatient(contains(string(tempPatient),'LatBHigh')) = 'LatBHigh';
tempPatient(contains(string(tempPatient),"0'HC")) = 'HC';
tempPatient(contains(string(tempPatient),"'HD - ")) = 'HD';
idx = double.empty(length(tempPatient),0);
idx(tempPatient == 'HC') = 1;
idx(tempPatient == 'HD') = 2;
idx(tempPatient == 'LatBLow') = 3;
idx(tempPatient == 'LatBHigh') = 4;
% Perform PCA on the data
[coeff,score,latent] = pca(z);
Y = score(:,1:2);

figure
gscatter(Y(:,1), Y(:,2),idx);
% hold on
% plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
legend('HC','HD','LatBLow','LatBHigh');
title('PCA Clustering');

Y = tsne(z, 'NumDimensions', num_dims);
figure
gscatter(Y(:,1), Y(:,2),idx);
% hold on
% plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
legend('HC','HD','LatBLow','LatBHigh');
title('TSNE Clustering');
%% 2d look at features by Well
close all
clc
num_dims = 2; % Reduce the data to 2 dimensions
% labels = regexp(string(uniquePatientsWell),'\d*','Match', 'once');
uniquePatientsWell(contains(string(uniquePatientsWell),"'HDMedia")) = 'HDMedia';
uniquePatientsWell(contains(string(uniquePatientsWell),"0'HCMedia")) = 'HCMedia';
uniquePatientsWell(contains(string(uniquePatientsWell),'LatB')) = 'LatB';
uniquePatientsWell(contains(string(uniquePatientsWell),"0'HC")) = 'HC';
uniquePatientsWell(contains(string(uniquePatientsWell),"'HD - ")) = 'HD';
idx = double.empty(length(uniquePatientsWell),0);
idx(uniquePatientsWell == 'HC') = 1;
idx(uniquePatientsWell == 'HD') = 2;
idx(uniquePatientsWell == 'LatB') = 3;
% Perform PCA on the data
[coeff,score,latent] = pca(zWell);
Y = score(:,1:2);
matrix = cell2mat(newcolorsRGB(1:4));
matrix = reshape(matrix, 3, []).';
figure;
hold on;
for i = 1:max(idx)
    scatter(Y(idx==i,1), Y(idx==i,2), 'filled', 'MarkerFaceColor', matrix(i, :), 'MarkerEdgeColor', matrix(i, :), 'SizeData', 50);
end
hold off;
% Get unique group indices
groupIndices = unique(idx);
% Add mean for each group
hold on
for i = 1:length(groupIndices)
    groupIndex = groupIndices(i);
    groupData = Y(idx == groupIndex, :);
    groupMean = mean(groupData);
    
    % Add a star for the mean of the group
% Plot the outline
plot(groupMean(1), groupMean(2), ...
    'MarkerEdgeColor', 'k', 'LineWidth', 2.5, 'Marker', 'p', 'MarkerSize', 11); % Slightly bigger size

hold on

% Plot the fill
plot(groupMean(1), groupMean(2), ...
    'MarkerEdgeColor', matrix(i, :), 'LineWidth', 1.5, 'Marker', 'p', 'MarkerSize', 10); % Original size
end
hold off
% hold on
% plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
legend('HC','HD','LatB  60nM', 'location', 'southwest');
title('PCA Clustering');
grid on;
    axis square;
    box on
    g = gcf;
g.WindowState = 'maximized';
saveas(g, strcat(folderCurrent, '\PCA2d.png'))

% Y = tsne(zWell, 'NumDimensions', num_dims);
% figure
% gscatter(Y(:,1), Y(:,2),idx);
% hold on
% plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
% legend('HC','HD','LatBLow','LatBHigh');
% title('TSNE Clustering');
%% 2d look at features
close all
clc
% labels = regexp(string(uniquePatientsWell),'\d*','Match', 'once');
uniquePatientsWell(contains(string(uniquePatientsWell),"'HDMedia")) = 'HDMedia';
uniquePatientsWell(contains(string(uniquePatientsWell),"0'HCMedia")) = 'HCMedia';
uniquePatientsWell(contains(string(uniquePatientsWell),'LatBLow')) = 'LatBLow';
uniquePatientsWell(contains(string(uniquePatientsWell),'LatBHigh')) = 'LatBHigh';
uniquePatientsWell(contains(string(uniquePatientsWell),"0'HC")) = 'HC';
uniquePatientsWell(contains(string(uniquePatientsWell),"'HD - ")) = 'HD';


% Perform PCA on the data
[coeff,score,latent] = pca(zWell);

% Keep the first two principal components
pc1 = score(:,1);
pc2 = score(:,2);

k = 4; 
Y = score(:,1:2);

% perplexity = 30; % Replace 30 with the desired perplexity value
% num_dims = 2; % Reduce the data to 2 dimensions
% 
% % Perform t-SNE on the data
% Y = tsne(z, 'NumDimensions', num_dims, 'Perplexity', perplexity);
% 

% Cluster the data using k-means
[idx,C] = kmeans(Y,k);
for ki=1:1:k
    disp('KMEANS');
    disp(strcat('HC members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HC")))));
    disp(strcat('HD members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HD")))));
    disp(strcat('LatB High members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "LatBHigh")))));
    disp(strcat('LatB Low members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "LatBLow")))));
    disp(strcat('HDMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HDMedia")))));
    disp(strcat('HCMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HCMedia")))));
end
disp('***************************');

% Plot the clustered data
figure
gscatter(Y(:,1), Y(:,2),idx);
hold on
plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Centroids');
title('KMeans Clustering');
%%

eps = 0.5; % Replace 0.5 with the desired epsilon value
MinPts = 5; % Replace 5 with the desired minimum number of points in a cluster

% Cluster the data using DBSCAN
[idx, C] = dbscan(Y, eps, MinPts);
for ki=1:1:length(unique(idx))
    disp('DBSCAN Clustering');
    disp(strcat('HC members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HC")))));
    disp(strcat('HD members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HD")))));
    disp(strcat('MitoQ members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "MitoQ")))));
    disp(strcat('HDMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HDMedia")))));
    disp(strcat('HCMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HCMedia")))));
end
disp('***************************');

% Plot the clustered data
figure
gscatter(Y(:,1), Y(:,2),idx);
title(['DBSCAN Clustering (\epsilon = ' num2str(eps) ', MinPts = ' num2str(MinPts) ')']);
%
Z = linkage(Y, 'ward');
% dendrogram(Z);

% Choose the number of clusters and assign them
idx = cluster(Z, 'Maxclust', k);
for ki=1:1:k
    disp('Hierarchical Clustering');
    disp(strcat('HC members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HC")))));
    disp(strcat('HD members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HD")))));
    disp(strcat('MitoQ members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "MitoQ")))));
    disp(strcat('HDMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HDMedia")))));
    disp(strcat('HCMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HCMedia")))));
end
disp('***************************');

% Plot the clustered data
figure
gscatter(Y(:,1), Y(:,2),idx);
title(['Hierarchical Clustering (k = ' num2str(k) ')']);

idx = spectralcluster(Y, k, 'Distance', 'euclidean');
for ki=1:1:k
    disp('Spectral Clustering');
    disp(strcat('HC members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HC")))));
    disp(strcat('HD members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HD")))));
    disp(strcat('MitoQ members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "MitoQ")))));
    disp(strcat('HDMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HDMedia")))));
    disp(strcat('HCMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HCMedia")))));
end
disp('***************************');

% Plot the clustered data
figure
gscatter(Y(:,1), Y(:,2),idx);
title(['Spectral Clustering (k = ' num2str(k) ')']);

Z = linkage(Y,'ward');
idx = cluster(Z,'Maxclust',k);
for ki=1:1:k
    disp('Agglomerative Clustering');
    disp(strcat('HC members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HC")))));
    disp(strcat('HD members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HD")))));
    disp(strcat('MitoQ members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "MitoQ")))));
    disp(strcat('HDMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HDMedia")))));
    disp(strcat('HCMedia members in cluster ', {' '}, string(ki), {' '}, 'is ', {' '}, string(sum(ismember(uniquePatientsWell(idx == ki), "HCMedia")))));
end
disp('***************************');

% Plot the clustered data
figure
gscatter(Y(:,1), Y(:,2),idx);
title(['Agglomerative Clustering (k = ' num2str(k) ')']);
%% Clusters Wells
runInsiders = false;
if runInsiders
    folderCurrent = folderCurrent + 'Insiders\';
    mkdir(folderCurrent);
end

% temp_CAPPatient = string(data.CAPPatient(combinedIndex));

cMap = [linspace(0,1,128)'*[1,1], ones(128,1)]; % The blue range
cMap = [cMap; rot90(cMap,2)];                   % add the red range
cMap = (cMap);
    
rowLabels = string(temp_cellType);
for i=1:size(rowLabels, 1)
    rowLabels(i) = rowLabels(i) + '_' + string(i) + '_' + string(temp_cellType(i));
end
cg = clustergram(normalized_features, ...
                             'ColumnLabels', feature_names,...
                             'RowLabels', rowLabels,...
                             'RowPdist', 'cityblock',...
                             'ColumnPdist', 'cityblock',...
                             'Linkage', 'ward',...
                             'Cluster', 'all',...
                             'Dendrogram', 75,'Colormap', cMap,...
                           "Standardize", 'column');
%% Features
newcolors = {'#4285F4','#0F9D58','#F4B400','#DB4437' '#330072', '#016773'};
newcolorsRGB = {[66 133 244]/255,[15 157 88]/255,[244 180 0]/255,[219 68 55]/255, [51 0 114]/255};
folderCurrent = "C:\Cell_Migration_MATLAB\NatureAviv\Compunds2\Features\LatB\Plots\";
mkdir(folderCurrent);

   clusterNames = string(unique(data.CellType(combinedIndex)));

if ~isempty(clusterNames)

    vNumbersAll = [];
    normedAll = [];
    vMeansAll = [];
    importanceAll = [];
    importanceNormed = [];

zClusters = z;
normedzClusters = rescale(zClusters, 'InputMin',min(zClusters),'InputMax',max(zClusters));

    for iBig=1:1:length(clusterNames)
    starts = iBig;
    ends = iBig+0;
    if ends > length(clusterNames)
        ends = length(clusterNames); 
    end
for myGroup=starts:ends

        currentGroupName = clusterNames{myGroup};
        if strfind(currentGroupName, 'HC') > 0
            myColorO = newcolors{1};
            myTypeOuter = 'HC';
        elseif strfind(currentGroupName, 'HD') > 0
            myColorO = newcolors{3};
            myTypeOuter = 'HD';
        elseif strfind(currentGroupName, 'LatBLow') > 0
            myColorO = newcolors{4};
            myTypeOuter = 'LatBLow';            
        elseif strfind(currentGroupName, 'LatBHigh') > 0
            myColorO = newcolors{5};
            myTypeOuter = 'LatBHigh';    
        end
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(temp_cellType == clusterNames(myGroup)) = 1;
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
    temp_model = fitcensemble(z,aNumbers, 'Method','bag','NumLearningCycles', 100);
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
        for lookedGroup=1:length(clusterNames)
        number = temp_cellType == clusterNames(myGroup);
        aFeature = [aFeature; z(number, featureToCheck)];
%             aPol = [aPol; xys(number)];
%             vColors = repmat({myColor},size(xys(number),1),1);
%             aColors = [aColors; vColors];
%             vType = repmat({myType},size(xys(number),1),1);
%             aType = [aType; vType];
        end
        nexttile;
        [f1,xi1] = ksdensity(aFeature, 'Kernel','normal', 'Support',[floor(min(aFeature)) ceil(max(aFeature))], 'NumPoints',size(aFeature,1));
        area(xi1,f1, 'FaceColor', myColorO, 'FaceAlpha', 0.5)
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
%     normedZ = rescale(z, 'InputMin',min(z),'InputMax',max(z));
%     r = figure;
%     normed = normedZ(vNumbers, :);
% %     normed = rescale(normed, 'InputMin',min(z),'InputMax',max(z));
%     normed = mean(normed);
%     zerosS = 0*ones(1,size(normed,2));
%     onesS = ones(1,size(normed,2));
%     L = [zerosS;onesS];
%     spider_plot(normed,'AxesLimits',L, ...
%         'AxesPrecision', onesS,...
%         'FillOption', {'on'},...
%         'Color', myColorO, ...
%         'FillTransparency', 0.15,...
%         'AxesLabels', feature_names,...
%         'AxesLabelsColors', ismember(string(feature_names), string(feature_importance_top)), ...
%         'AxesInterval', 5);
%     title(strcat('Spider Plot (by All Data) for',{' '}, string(clusterNames{myGroup})));
%     exportgraphics(r,strcat(folderCurrent, 'g_', string(myGroup),'_RadarTotal.jpg'),'Resolution',500)
%     
%     myText = sprintf("I have a cluster with %f Dp, %f Dtot, %f Pp, %f Psi, %f Pnp, %f Dnp, %f MSD, %f Sp, %f Snp. What could it mean?", normed(1), normed(2), normed(3), normed(4), normed(5), normed(6), normed(7), normed(8), normed(9));
%     disp(myText);
    %
    % Radar Local
    r = figure;
    vNumbers = temp_cellType == clusterNames(myGroup);
    normed = z(vNumbers, :);
    normed = mean(normed);
    normed = rescale(normed, 'InputMin',-1,'InputMax',1);
%     normedAll = [normedAll; normed];
    zerosS = 0*ones(1,size(normed,2));
    onesS = ones(1,size(normed,2));
    L = [zerosS;onesS];
    spider_plot(normed,'AxesLimits',L, ...
        'AxesPrecision', onesS, ...
        'FillOption', {'on'},...
        'Color', myColorO, ...
        'FillTransparency', 0.15,...
        'AxesLabels', feature_names,...
        'AxesLabelsColors', ismember(string(feature_names), string(feature_importance_top)), ...
        'AxesInterval', 5);
    title(strcat('Spider Plot (by Clusters Data) for',{' '}, string(clusterNames{myGroup})));
    exportgraphics(r,strcat(folderCurrent, 'h_', string(myGroup),'_RadarClusters.jpg'),'Resolution',500)
    
    myText = sprintf("I have a cluster with %f Dp, %f Dtot, %f Pp, %f Psi, %f Pnp, %f Dnp, %f MSD, %f Sp, %f Snp. What could it mean?", normed(1), normed(2), normed(3), normed(4), normed(5), normed(6), normed(7), normed(8), normed(9));
    disp(myText);
end
    end
end
                       %%

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
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-HC")};
                myCellHC = [myCellHC; strcat('Cluster-',string(groupIndexed), "-HC"), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellHCCAPS = [myCellHCCAPS; strcat('Cluster-',string(groupIndexed), "-HC"), CAPScoresString];

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
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-Severe")};
                myCellS = [myCellS; strcat('Cluster-',string(groupIndexed), "-Severe"), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellSCAPS = [myCellSCAPS; strcat('Cluster-',string(groupIndexed), "-Severe"), CAPScoresString];
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
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-Mild")};
                myCellM = [myCellM; strcat('Cluster-',string(groupIndexed), "-Mild"), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellMCAPS = [myCellMCAPS; strcat('Cluster-',string(groupIndexed), "-Mild"), CAPScoresString];

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
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-Premanifest")};
                myCellR = [myCellR; strcat('Cluster-',string(groupIndexed), "-Premanifest"), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellRCAPS = [myCellRCAPS; strcat('Cluster-',string(groupIndexed), "-Premanifest"), CAPScoresString];
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
                clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed), "-HGPS")};
                myCellP = [myCellP; strcat('Cluster-',string(groupIndexed), "-HGPS"), countHC, countRisk, countMild, countSevere, countProgeria];
                myCellPCAPS = [myCellPCAPS; strcat('Cluster-',string(groupIndexed), "-HGPS"), CAPScoresString];
            end
        end
        
        disp(strcat("Total_Cells_Is_", string(countMild + countRisk + countSevere + countHC + countProgeria)));

        end
end
    
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

%
if ~isempty(clusterLabels)
rm = struct('GroupNumber',clusterLabels,'Annotation',clusterNames,...
    'Color',CMCell, 'FontSize', 5);
set(cg,'RowGroupMarker',rm)
cgf = plot(cg);
set(cgf,'FontSize',14, 'fontweight','bold')
screensize = get( groot, 'Screensize' );
set(cgf, 'DefaultFigurePosition', screensize);
saveas(cgf,strcat(folderCurrent, '\c_Clustergram.jpg'))
end
%%
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
    saveas(g, strcat(folderCurrent, '\e_Bar_' + string(item) +'.png'))
end

%%
close all force
%$%$ Draw
       clusterNames = string(unique(data.CellType(combinedIndex)));

    for iBig=1:1:4
    starts = iBig;
    ends = iBig+0;
    if ends > 4
        ends = 4; 
    end
    fig = figure('Position', get(0, 'Screensize'));
    t = tiledlayout(3,1, 'TileSpacing','Compact');
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
            plot(xy(:,1)+counter_x*space,xy(:,2)+counter_y*space,'Color',col, 'LineWidth', 2);
            counter_x = counter_x + 1;
            if mod(counter_x,10)==0
                counter_x = 0;
                counter_y = counter_y + 1;
            end
        end
        set(gca,'xtick',[],'ytick',[]);
        title(strcat("Trajectories for ", string(clusterNames{myGroup})));
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
        title(strcat("Sunplot Trajectories for ", string(clusterNames{myGroup})));
        axis square
        hold off;
        set(gca,'xtick',[],'ytick',[]);
        ylim([-350 350])
        xlim([-350 350])
        
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
        plot(xy1,xy2,'-','color',myColorO, 'linewidth',2);
        title(strcat("Merged Trajectories for ", string(clusterNames{myGroup})));
        axis square
        hold off;
        set(gca,'xtick',[],'ytick',[]);
        ylim([-30 30])
        xlim([-30 30])

    end
    exportgraphics(t,strcat(folderCurrent, 'e_', string(myGroup),'_TRJ.jpg'),'Resolution',500)
    end
close all force

%% Features
   clusterNames = string(unique(data.CellType(combinedIndex)));
    vNumbersAll = [];
    normedAll = [];
    vMeansAll = [];
    importanceAll = [];
    importanceNormed = [];

zClusters = z;
normedzClusters = rescale(zClusters, 'InputMin',min(zClusters),'InputMax',max(zClusters));

    for iBig=1:1:length(clusterNames)
    starts = iBig;
    ends = iBig+0;
    if ends > length(clusterNames)
        ends = length(clusterNames); 
    end
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
        if strfind(currentGroupName, 'HC') > 0
            myColorO = newcolors{1};
            myTypeOuter = 'HC';
        elseif strfind(currentGroupName, 'HD') > 0
            myColorO = newcolors{3};
            myTypeOuter = 'HD';
        elseif strfind(currentGroupName, 'LatBLow') > 0
            myColorO = newcolors{4};
            myTypeOuter = 'LatBLow';
        elseif strfind(currentGroupName, 'LatBHigh') > 0
            myColorO = newcolors{5};
            myTypeOuter = 'LatBHigh';
        end
        aNumbers = zeros(size(z, 1), 1);
        aNumbers(temp_cellType == clusterNames(myGroup)) = 1;
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
        for lookedGroup=1:length(clusterNames)
        number = temp_cellType == clusterNames(myGroup);
        aFeature = [aFeature; z(number, featureToCheck)];
%             aPol = [aPol; xys(number)];
%             vColors = repmat({myColor},size(xys(number),1),1);
%             aColors = [aColors; vColors];
%             vType = repmat({myType},size(xys(number),1),1);
%             aType = [aType; vType];
        end
        nexttile;
        [f1,xi1] = ksdensity(aFeature, 'Kernel','normal', 'Support',[floor(min(aFeature)) ceil(max(aFeature))], 'NumPoints',size(aFeature,1));
        area(xi1,f1, 'FaceColor', myColorO, 'FaceAlpha', 0.5)
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
%     normedZ = rescale(z, 'InputMin',min(z),'InputMax',max(z));
%     r = figure;
%     normed = normedZ(vNumbers, :);
% %     normed = rescale(normed, 'InputMin',min(z),'InputMax',max(z));
%     normed = mean(normed);
%     zerosS = 0*ones(1,size(normed,2));
%     onesS = ones(1,size(normed,2));
%     L = [zerosS;onesS];
%     spider_plot(normed,'AxesLimits',L, ...
%         'AxesPrecision', onesS,...
%         'FillOption', {'on'},...
%         'Color', myColorO, ...
%         'FillTransparency', 0.15,...
%         'AxesLabels', feature_names,...
%         'AxesLabelsColors', ismember(string(feature_names), string(feature_importance_top)), ...
%         'AxesInterval', 5);
%     title(strcat('Spider Plot (by All Data) for',{' '}, string(clusterNames{myGroup})));
%     exportgraphics(r,strcat(folderCurrent, 'g_', string(myGroup),'_RadarTotal.jpg'),'Resolution',500)
%     
%     myText = sprintf("I have a cluster with %f Dp, %f Dtot, %f Pp, %f Psi, %f Pnp, %f Dnp, %f MSD, %f Sp, %f Snp. What could it mean?", normed(1), normed(2), normed(3), normed(4), normed(5), normed(6), normed(7), normed(8), normed(9));
%     disp(myText);
    %
    % Radar Local
    r = figure;
    vNumbers = temp_cellType == clusterNames(myGroup);
    normed = z(vNumbers, :);
    normed = mean(normed);
    normed = rescale(normed, 'InputMin',min(z),'InputMax',max(z));
%     normedAll = [normedAll; normed];
    zerosS = 0*ones(1,size(normed,2));
    onesS = ones(1,size(normed,2));
    L = [zerosS;onesS];
    spider_plot(normed,'AxesLimits',L, ...
        'AxesPrecision', onesS, ...
        'FillOption', {'on'},...
        'Color', myColorO, ...
        'FillTransparency', 0.15,...
        'AxesLabels', feature_names,...
        'AxesLabelsColors', ismember(string(feature_names), string(feature_importance_top)), ...
        'AxesInterval', 5);
    title(strcat('Spider Plot (by Clusters Data) for',{' '}, string(clusterNames{myGroup})));
    exportgraphics(r,strcat(folderCurrent, 'h_', string(myGroup),'_RadarClusters.jpg'),'Resolution',500)
    
    myText = sprintf("I have a cluster with %f Dp, %f Dtot, %f Pp, %f Psi, %f Pnp, %f Dnp, %f MSD, %f Sp, %f Snp. What could it mean?", normed(1), normed(2), normed(3), normed(4), normed(5), normed(6), normed(7), normed(8), normed(9));
    disp(myText);
end
    end
    %%
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

%%
% Prepare 1-dimensional displacement trajectories
xys = xysMain;

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
scatter(xAxis, one_dim_xys(r,1:74),sz, c, 'filled'); hold on;
colorbar
colormap jet
plot(xAxis, one_dim_xys(r,1:74),'-','color','k', 'linewidth',0.5); hold off;
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

temp_cellType = data.CellType;
% temp_CAPPatient = string(data.CAPPatient(idxComb));

rowLabels = string(temp_cellType);
for i=1:size(rowLabels, 1)
    rowLabels(i) = rowLabels(i) + '_' + string(i) + '_' + string(temp_cellType(i));
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
GI_One_Dim = [g1; g2; g3; g5];
GINums = [4615 4603 4570 4614];
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
% uniqueCAPS = natsort(unique(regexp(string(data.CAPPatient2(combinedIndex)),'\d*','Match', 'once')));
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
    countSevere = sum(contains(string(currentGroupInfo.RowNodeNames), "HD"));
    countMild = sum(contains(string(currentGroupInfo.RowNodeNames), "MitoQ"));
    countRisk = sum(contains(string(currentGroupInfo.RowNodeNames), "Premanifest"));
    countHC = sum(contains(string(currentGroupInfo.RowNodeNames), "HC"));
    countProgeria = sum(contains(string(currentGroupInfo.RowNodeNames), "HGPS"));
    
%     CAPScores = containers.Map;
%     for uC=1:length(uniqueCAPS)
%         CAPScores(uniqueCAPS(uC)) = 0;
%     end
    for item=1:length(currentGroupInfo.RowNodeNames)
        currString = string(currentGroupInfo.RowNodeNames(item));
        split_string = strsplit(currString, '_');
        resultCAP = regexp(split_string{3},'\d*','Match', 'once');
%         CAPScores(resultCAP) = CAPScores(resultCAP) + 1;
    end
%     [~,sortedCAPScoresIndexes] = natsort(CAPScores.keys);
%     CAPScoresString = string(CAPScores.values);
%     CAPScoresString = CAPScoresString(sortedCAPScoresIndexes);
    
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
%                 myCellHCCAPS = [myCellHCCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];

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
%                 myCellSCAPS = [myCellSCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];
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
%                 myCellMCAPS = [myCellMCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];

           end
        end
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
%                 clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed))};
%                 myCellR = [myCellR; strcat('Cluster-',string(groupIndexed)), countHC, countRisk, countMild, countSevere, countProgeria];
%                 myCellRCAPS = [myCellRCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];
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
%                 clusterNames(end+1) = {strcat('Cluster-',string(groupIndexed))};
%                 myCellP = [myCellP; strcat('Cluster-',string(groupIndexed)), countHC, countRisk, countMild, countSevere, countProgeria];
%                 myCellPCAPS = [myCellPCAPS; strcat('Cluster-',string(groupIndexed)), CAPScoresString];
%             end
%         end
        
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
rm = struct('GroupNumber',clusterLabels,'Annotation',clusterNames,...
    'Color',CMCell, 'FontSize', 5);
set(cg,'RowGroupMarker',rm)
cgf = plot(cg);
set(cgf,'FontSize',14, 'fontweight','bold')
screensize = get( groot, 'Screensize' );
set(cgf, 'DefaultFigurePosition', screensize);
% saveas(cgf,strcat(folderCurrent, '\c_Clustergram.jpg'))
end
%%
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
   clusterNames = string(unique(data.CellType(combinedIndex)));
    for iBig=1:1:length(clusterNames)
    starts = iBig;
    ends = iBig+0;
    if ends > 4
        ends = 4;
    end
%     fig = figure('Position', get(0, 'Screensize'));
%     t = tiledlayout(3,1, 'TileSpacing','Compact');
    for myGroup=starts:ends
        aPol = [];
        aColors = [];
        aType = [];
        RiskCount = 0;
        MildCount = 0;
        SevereCount = 0;
        HCCount = 0;
        ProgeriaCount = 0;
        number = data.CellType(combinedIndex) == clusterNames(myGroup);
        aPol = [aPol one_dim_xys(number,:)'];
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
            title(['Activity Profile for ' clusterNames(myGroup)])
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
    saveas(g, strcat(folderCurrent, '\Ag_Activity_' + string(iBig) + '.png'))

    end
%% Lags and Trains
% binarize
close all
meanArray = mean(one_dim_xys,2);
stdArray = std(one_dim_xys,0,2);
binarized = one_dim_xys>meanArray+stdArray;
TrainSum = [];
LagsSum = [];
TrainCV = [];
LagsCV = [];
grp1 = [];
labels = [];

   clusterNames = string(unique(data.CellType(combinedIndex)));
    for iBig=1:1:length(clusterNames)
    starts = iBig;
    ends = iBig+0;
    if ends > 4
        ends = 4;
    end
    for myGroup=starts:ends
        labels = clusterNames;
        TrainsPol = [];
        LagsPol = [];
        number = data.CellType(combinedIndex) == clusterNames(myGroup);
        TrainsPol = [TrainsPol; binarized(number,:)];
        LagsPol = [LagsPol; 1-binarized(number,:)];
        TrainSum = [TrainSum, sum(TrainsPol, 2)'];
        LagsSum = [LagsSum, sum(LagsPol, 2)'];
        TrainCV = [TrainCV, (std(TrainsPol,0, 2)./mean(TrainsPol, 2))'];
        LagsCV = [LagsCV, (std(LagsPol,0, 2)./mean(LagsPol, 2))'];
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
    g = gcf;
    saveas(g, strcat(folderCurrent, '\Ah_Trains_Length.png'))
    figure()
    boxplot(LagsSum, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', labels)
    title('Boxplot for Lags');
    ylabel('Lags Length');
    ylim([0 max(LagsSum)*1.05]);
    g = gcf;
    saveas(g, strcat(folderCurrent, '\Ah_Lags_Length.png'))
    figure()
    boxplot(TrainCV, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', labels)
    title('Boxplot for CV Trains');
    ylabel('CV Train Length');
    ylim([0 max(TrainCV)*1.05]);
    g = gcf;
    saveas(g, strcat(folderCurrent, '\Ah_CVTrains_Length.png'))
    figure()
    boxplot(LagsCV, grp1, 'colorgroup', grp1, 'boxstyle', 'outline', 'labels', labels)
    title('Boxplot for CV Lags');
    ylabel('CV Lags Length');
    ylim([0 max(LagsCV)*1.05]);
    g = gcf;
    saveas(g, strcat(folderCurrent, '\Ah_CVLags_Length.png'))

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
        M(iA, iM) = num_common;% / length(currentClusterActivity);
    end
end
heatmap(string(clusterNamesMotility),string(clusterNamesActivity), M/max(max(M)), 'Colormap', parula);
title('Abundance of Activity cluster per Motility cluster')
g = gcf;
saveas(g, strcat(folderCurrent, '\Aj_heatmap_motility.png'))
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
