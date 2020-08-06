
%import data from excel file (simple version)

fileName='test import.xlsx';
[~,sheetName] = xlsfinfo(fileName);
behavData=cell(numel(sheetName),1);
for dayNum=1:numel(sheetName)
    [behavData{dayNum} strings] = xlsread(fileName,sheetName{dayNum}); 
end



