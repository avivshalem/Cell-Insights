VarName1All = [VarName1];
VarName2All = [VarName2];
VarName3All = [VarName3];
VarName4All = [VarName4];
VarName5All = [VarName5];
VarName6All = [VarName6];
VarName7All = [VarName7];

%% Find when the microscope stopped imaging
min_frame = 75;
% end_idx = find(FrameNumber == min_frame);
% end_idx = find(VarName6 == min_frame);
end_idx = find(VarName6All == min_frame);
new_CELLTYPE = [];
new_CELLID = [];
new_VarName3 = [];
new_ycoordinatepixel = [];
new_xcoordinatepixel = [];
new_CellID = [];
new_Distancefromoriginpixel = [];
new_Distancefromlastpointpixel = [];
new_Instantaneousspeedpixelminute = [];
new_Anglefromorigindegree = [];
new_Anglefromlastpointdegree = [];
new_Framenumber = [];

for i=1:length(end_idx)
    temp_idx = end_idx(i);
    % find temp data
%     temp_Framenumber = FrameNumber(temp_idx-min_frame+1:temp_idx);
%     temp_CELLTYPE = CELLTYPE(temp_idx-min_frame+1:temp_idx);
%     temp_CELLID = CELLID(temp_idx-min_frame+1:temp_idx);
% %     temp_VarName3 = VarName3(temp_idx-min_frame+1:temp_idx);
%     temp_ycoordinatepixel = ycoordinatepixel(temp_idx-min_frame+1:temp_idx);
%     temp_xcoordinatepixel = xcoordinatepixel(temp_idx-min_frame+1:temp_idx);
%     temp_CellID = CellID(temp_idx-min_frame+1:temp_idx);
%     %temp_Distancefromoriginpixel = Distancefromoriginpixel(temp_idx-min_frame+1:temp_idx);
%     %temp_Distancefromlastpointpixel = Distancefromlastpointpixel(temp_idx-min_frame+1:temp_idx);
%     %temp_Instantaneousspeedpixelminute = Instantaneousspeedpixelminute(temp_idx-min_frame+1:temp_idx);
%     %temp_Anglefromorigindegree = Anglefromorigindegree(temp_idx-min_frame+1:temp_idx);
%     %temp_Anglefromlastpointdegree = Anglefromlastpointdegree(temp_idx-min_frame+1:temp_idx);
    
%     temp_Framenumber = VarName6(temp_idx-min_frame+1:temp_idx);
%     temp_CELLTYPE = VarName2(temp_idx-min_frame+1:temp_idx);
%     temp_CELLID = VarName3(temp_idx-min_frame+1:temp_idx);
    temp_Framenumber = VarName6All(temp_idx-min_frame+1:temp_idx);
    temp_CELLTYPE = VarName2All(temp_idx-min_frame+1:temp_idx);
    temp_CELLID = VarName3All(temp_idx-min_frame+1:temp_idx);
    %temp_VarName3 = Column1(temp_idx-min_frame+1:temp_idx);
%     temp_ycoordinatepixel = VarName5(temp_idx-min_frame+1:temp_idx);
%     temp_xcoordinatepixel = VarName4(temp_idx-min_frame+1:temp_idx);
%     temp_CellID = VarName7(temp_idx-min_frame+1:temp_idx);
    temp_ycoordinatepixel = VarName5All(temp_idx-min_frame+1:temp_idx);
    temp_xcoordinatepixel = VarName4All(temp_idx-min_frame+1:temp_idx);
    temp_CellID = VarName7All(temp_idx-min_frame+1:temp_idx);
    
    %insert data to new vectors
    new_Framenumber = [new_Framenumber; temp_Framenumber];
    new_CELLTYPE = [new_CELLTYPE; temp_CELLTYPE];
    new_CELLID = [new_CELLID; temp_CELLID];
%     new_VarName3 = [new_VarName3; temp_VarName3];
    new_ycoordinatepixel = [new_ycoordinatepixel; temp_ycoordinatepixel];
    new_xcoordinatepixel = [new_xcoordinatepixel; temp_xcoordinatepixel];
    new_CellID = [new_CellID; temp_CellID];
    %new_Distancefromoriginpixel = [new_Distancefromoriginpixel; temp_Distancefromoriginpixel];
    %new_Distancefromlastpointpixel = [new_Distancefromlastpointpixel; temp_Distancefromlastpointpixel];
    %new_Instantaneousspeedpixelminute = [new_Instantaneousspeedpixelminute; temp_Instantaneousspeedpixelminute];
    %new_Anglefromorigindegree = [new_Anglefromorigindegree; temp_Anglefromorigindegree];
    %new_Anglefromlastpointdegree = [new_Anglefromlastpointdegree; temp_Anglefromlastpointdegree];

end
%% Change Cell Id
flag = 1;
for i=1:min_frame:numel(new_CellID)
    new_CellID(i:i + min_frame - 1) = flag;
    %new_ycoordinatepixel(i:i + min_frame - 1) = new_ycoordinatepixel(i:i + min_frame - 1) - new_ycoordinatepixel(i);
    %new_xcoordinatepixel(i:i + min_frame - 1) = new_xcoordinatepixel(i:i + min_frame - 1) - new_xcoordinatepixel(i);
    flag = flag + 1;
end
%% Write to xlsx
T = table(new_CELLID, new_ycoordinatepixel, new_xcoordinatepixel,...
    new_Framenumber, new_CellID);
filename = 'C:\CellInsights\Tracking\Data_All_Information.xlsx';
writetable(T,filename,'Sheet',1)

T = table(new_CellID, new_Framenumber, new_xcoordinatepixel, new_ycoordinatepixel);
filename = 'C:\CellInsights\Tracking\Data_CellID_Frame_X_Y.xlsx';
writetable(T,filename,'Sheet',1)
 
%% MEAN SQUARED DISPLACEMENT 10
msd = get_MSD;

%% FITTING OF CELL TRAJECTORIES TO THE APRW MODEL 10
% [P, S, Positioning Error, R Squared. RMSE]
fit_APRW;

%% Fix ID's
CellID_unique = unique(new_CellID);
CELLTYPE_unique = categorical(zeros(size(CellID_unique)));
CELLID_unique = categorical(zeros(size(CellID_unique)));
%VarName3_unique = categorical(zeros(size(CellID_unique)));

for i=1:numel(CellID_unique)
    temp_CellID = CellID_unique(i);
    CELLTYPE_unique(i) = new_CELLTYPE(find(new_CellID == temp_CellID,1,'first'));
    CELLID_unique(i) = new_CELLID(find(new_CellID == temp_CellID,1,'first'));
    %VarName3_unique(i) = new_VarName3(find(new_CellID == temp_CellID,1,'first'));
end

%% Below are Extras, not needed
% 
% 
% 
% 
% 
%
%
%
%
%
%% Write to Mat
T = [new_CellID, new_Framenumber, new_xcoordinatepixel, new_ycoordinatepixel];

filename = 'Data.mat';
save(filename, 'T')

%%
VarName1All = [VarName1; VarName8];
VarName2All = [VarName2; VarName9];
VarName3All = [VarName3; VarName10];
VarName4All = [VarName4; VarName11];
VarName5All = [VarName5; VarName12];
VarName6All = [VarName6; VarName13];
VarName7All = [VarName7; VarName14];
%%
VarName1All = [VarName1; VarName8; VarName15];
VarName2All = [VarName2; VarName9; VarName16];
VarName3All = [VarName3; VarName10; VarName17];
VarName4All = [VarName4; VarName11; VarName18];
VarName5All = [VarName5; VarName12; VarName19];
VarName6All = [VarName6; VarName13; VarName20];
VarName7All = [VarName7; VarName14; VarName21];
%%
% VarName1All = [VarName1; VarName8; VarName15; VarName22];
% VarName2All = [VarName2; VarName9; VarName16; VarName23];
% VarName3All = [VarName3; VarName10; VarName17; VarName24];
% VarName4All = [VarName4; VarName11; VarName18; VarName25];
% VarName5All = [VarName5; VarName12; VarName19; VarName26];
% VarName6All = [VarName6; VarName13; VarName20; VarName27];
% VarName7All = [VarName7; VarName14; VarName21; VarName28];

% VarName1All = [VarName1; VarName8; VarName15; VarName22; VarName29; VarName36];
% VarName2All = [VarName2; VarName9; VarName16; VarName23; VarName30; VarName37];
% VarName3All = [VarName3; VarName10; VarName17; VarName24; VarName31; VarName38];
% VarName4All = [VarName4; VarName11; VarName18; VarName25; VarName32; VarName39];
% VarName5All = [VarName5; VarName12; VarName19; VarName26; VarName33; VarName40];
% VarName6All = [VarName6; VarName13; VarName20; VarName27; VarName34; VarName41];
% VarName7All = [VarName7; VarName14; VarName21; VarName28; VarName35; VarName42];
VarName1All = [VarName1; VarName8; VarName15; VarName22; VarName29; VarName36; VarName43];
VarName2All = [VarName2; VarName9; VarName16; VarName23; VarName30; VarName37; VarName44];
VarName3All = [VarName3; VarName10; VarName17; VarName24; VarName31; VarName38; VarName45];
VarName4All = [VarName4; VarName11; VarName18; VarName25; VarName32; VarName39; VarName46];
VarName5All = [VarName5; VarName12; VarName19; VarName26; VarName33; VarName40; VarName47];
VarName6All = [VarName6; VarName13; VarName20; VarName27; VarName34; VarName41; VarName48];
VarName7All = [VarName7; VarName14; VarName21; VarName28; VarName35; VarName42; VarName49];