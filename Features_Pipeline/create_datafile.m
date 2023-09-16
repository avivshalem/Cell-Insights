% clear all; close all;
%% Loading Data file APRW
DataFileMitoQ = struct();

DataFileMitoQ.CellID = CellID;
DataFileMitoQ.Dnp = Dnp;
DataFileMitoQ.Dp = Dp;
DataFileMitoQ.Dtot = Dtot;
DataFileMitoQ.MSD10 = MSD10;
DataFileMitoQ.Pnp = Pnp;
DataFileMitoQ.Snp = Snp;
DataFileMitoQ.Pp = Pp;
DataFileMitoQ.Sp = Sp;
DataFileMitoQ.SP10 = SP10;
DataFileMitoQ.Psi = Psi;
DataFileMitoQ.CellType = CELLTYPE;
%DataFile20122021.VarName3 = VarName3;
DataFileMitoQ.Patient = PatientCELLID;
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
DataFileMitoQ.CAP = CAP;
DataFileMitoQ.CAPPatient = CAPPatient;
DataFileMitoQ.CAPPatient2 = CAPPatient2;



save DataFileMitoQ DataFileMitoQ

%% trajectories (IF RAN ONCE, NOT NEEDED) 10
xys_full_MitoQ = sim_APRW;
%%
xys_TRJs_MitoQ = get_trajfile; %(DATA.mat)
%%
filename = 'xys_full_MitoQ.mat';
save(filename, 'xys_full_MitoQ')


filename2 = 'xys_TRJs_MitoQ.mat';
save(filename2, 'xys_TRJs_MitoQ')

%%

%% Loading Data file PRW
DataFile = struct();

DataFile.CellID = CellID;
% DataFile.Dnp = Dnp;
DataFile.Dp = Dp;
% DataFile.Dtot = Dtot;
DataFile.MSD10 = MSD10;
% DataFile.Pnp = Pnp;
% DataFile.Snp = Snp;
DataFile.Pp = Pp;
DataFile.Sp = Sp;
DataFile.SP10 = SP10;
% DataFile.Psi = Psi;
DataFile.CellType = CELLTYPE;
%DataFile20122021.VarName3 = VarName3;
DataFile.Patient = PatientCELLID;
% DataFile.Sex = Sex;
% DataFile.Sample = Sample;
% DataFile.Allele1 = Allele1;
% DataFile.Allele2 = Allele2;
% DataFile03012022_HD_HC.Family = Family;
% DataFile.Onset = Onset;
% DataFile.Delta = Delta;
% DataFile.MSD60 = MSD60;
%DataFile.SP60 = SP60;
DataFile.CAP = CAP;

save DataFile DataFile

%% Loading Data file APRW
DataFile = struct();

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
%DataFile20122021.VarName3 = VarName3;
DataFile.Patient = Patient;
DataFile.Well = Well;
DataFile.PatNW = PatientNoWell;
% DataFile.Sex = Sex;
% DataFile.Sample = Sample;
% DataFile.Allele1 = Allele1;
% DataFile.Allele2 = Allele2;
% DataFile03012022_HD_HC.Family = Family;
% DataFile.Onset = Onset;
% DataFile.Delta = Delta;
% DataFile.MSD60 = MSD60;
%DataFile.SP60 = SP60;
DataFile.CAP = CAP;
DataFile.CAPPatient = CAPPatient;


save DataFile DataFile

