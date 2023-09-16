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
