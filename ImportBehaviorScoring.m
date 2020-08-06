function behavData = ImportBehaviorScoring(workbookFile, sheetName)
% works just as well with [numbers strings misc] = xlsread(workbookFile,sheetName); 

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    startRow = 2;
    endRow = 81;
end

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + startRow(1) + ":D" + endRow(1);

% Specify column names and types
opts.VariableNames = ["Run", "Trial", "IsRunning", "IsSwiping"];
opts.VariableTypes = ["double", "double", "string", "string"];
opts = setvaropts(opts, [3, 4], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [3, 4], "EmptyFieldRule", "auto");

% Import the data
behavData = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:length(startRow)
    opts.DataRange = "A" + startRow(idx) + ":D" + endRow(idx);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    behavData = [behavData; tb]; %#ok<AGROW>
end

end