clear all; close all;
%% default path
default_path = 'C:\CellInsights\Tracking';

%% Loading Data file APRW
fileName = 'DataFile.mat';

DataFile = struct();

DataFile.CellID = CellID;
DataFile.Dnp = Dnp;
DataFile.Dp = Dp;
DataFile.Dtot = Dtot;
DataFile.MSD10 = MSD10;
DataFile.Pnp = Pnp;
DataFile.Snp = Snp;
DataFile.Pp = Pp;
DataFile.Sp = Sp;
DataFile.SP10 = SP10;
DataFile.Psi = Psi;
DataFile.CellType = CELLTYPE;
DataFile.Patient = PatientCELLID;
DataFile.CAP = CAP;
DataFile.CAPPatient = CAPPatient;
DataFile.CAPPatient2 = CAPPatient2;


fullFilePath = fullfile(default_path, fileName);
save(fullFilePath, 'DataFile');



% DataFile20122021.VarName3 = VarName3;
% DataFileWellAgeSex.Sex = Sex;
% DataFileWellAgeSex.Age = Age;
% DataFile.Sample = Sample;
% DataFile.Allele1 = Allele1;
% DataFile.Allele2 = Allele2;
% DataFile03012022_HD_HC.Family = Family;
% DataFile.Onset = Onset;
% DataFile.Delta = Delta;
% DataFile.MSD60 = MSD60;
%DataFile.SP60 = SP60;
%% TRJs APRW (IF RAN ONCE, NOT NEEDED) 10
xys_full = sim_APRW; %(APRW model fit.xlsx)
%% TRJs Full
xys_TRJs = get_trajfile; %(Data_CellID_Frame_X_Y.xlsx)
%%
filename = strcat(default_path, '\xys_full.mat');
save(filename, 'xys_full')


filename2 = strcat(default_path, '\xys_TRJs.mat');
save(filename2, 'xys_TRJs')

