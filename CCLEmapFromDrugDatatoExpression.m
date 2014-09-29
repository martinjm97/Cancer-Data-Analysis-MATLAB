% Reads in the gene expression table.
expression = readtable('CCLE_Expression_Entrez_2012-09-29.gct',...
    'FileType','text',...       Helps to read in this file type
    'Delimiter','\t',...
    'HeaderLines',2,...         Removes extra information
    'ReadVariableNames',true,...Makes columns
    'ReadRowNames',true);       % Makes rows


%% Reads in the gene expression table.

drugData = readtable('CCLE_NP24.2009_Drug_data_2012.02.20.csv',...
    'ReadVariableNames', true,... Allows for manual naming of variables
    'ReadRowNames', false,...     Names rows
    'Delimiter',',',...         Tells how to read data
    'TreatAsEmpty','NA',...     
    'Format','%s%s%s%s%q%q%q%f%s%f%f%f%f');    % fixes string error in table

%%
% matlab.lang.makeValidName is called during the construction of the table,
% change the names in the CCLE field of the drug data to match the variable
% names in the expression data

drugData.CCLECellLineName = matlab.lang.makeValidName(drugData.CCLECellLineName);

%% Finds all of the drug data for Topotecan

tempDrugData = drugData(strcmp(drugData.Compound,'Topotecan'),:);
size(tempDrugData)

% 504 cell lines were exposed to Topotecan

%%

compoundNames = unique(drugData.Compound);

for iter = 1:length(compoundNames)
    thisName = compoundNames{iter};
    numSamples(iter) = sum(strcmp(drugData.Compound,thisName)); 
    % sums a logical to find the number of samples for a given compound
end


%% One to one mapping?

% Checks that the cell lines that are in the expression set
% only have one match.

numStrMatches = zeros(size(tempDrugData,1),1);
for iter = 1:size(tempDrugData,1)
    numStrMatches(iter) = sum(strcmp(tempDrugData.CCLECellLineName{iter},expression.Properties.VariableNames));
end

% if there's a unique map between the two, then numStrMatches is 1 
% When the cell lines are not in the expression data it is 0.

%% Mapping
% Data Monging
% Maps from expression data to drug data

expression2Drug = zeros(size(tempDrugData,1),1);
isValidDrug = false(size(tempDrugData,1),1);

for iter = 1:size(tempDrugData,1)
    % compares the names in tempDrugData with those in CCLECellLineName
    % and retrieves their indicies.
    index = find(strcmp(tempDrugData.CCLECellLineName{iter},...
        expression.Properties.VariableNames));
    if isempty(index)
        % cell line is not in expression data
        expression2Drug(iter) = NaN;
        isValidDrug(iter) = false;
    else
        expression2Drug(iter) = index;
        isValidDrug(iter) = true;
    end
end

expression2Drug(isnan(expression2Drug)) = []; % remove indices that don't map

%%
% trims variables for analysis
trimmedTempDrugData = tempDrugData(isValidDrug,:);
trimmedExpression   = expression(:,expression2Drug);


%% Get the correct metagenes


% Loads metafeatures from a saved file such that the array has the exact same
% dimensions as expression.


run('importDataCCLEimproved.m');

%% meta = MetaWeights'*trimmedExpression{:,:};

% Calculates the metagenes using weights. 

meta = W'*trimmedExpression{:,:};   %makes metagenes
metaTransposed = meta';

%% Graphs the IC50 data

% Uses two metagenes graphed against eachother to find the two groups.
% Then graphs the two parts separately.

response = log10(trimmedTempDrugData.IC50_uM_); 
%log often used to for clinical trends

meta1 = metaTransposed (:,1); % takes the first and second metagenes
meta2 = metaTransposed (:,2);



figure(); hold on               % allows for multiple graphs
spacing = linspace(5,9);        % x-coordinates
height = 3/2*(spacing-5)+5;    % f(x)-coordinates made to separate groups
plot (spacing,height);     
plot (meta1,meta2,'o');

figure(); hold on
cdfplot(response( meta1>6 & meta2<5.7)); % thresholds taken from graph
cdfplot(response( meta2>5.7));

%% Graphs of the Amax

response = trimmedTempDrugData.Amax;

figure(); hold on
cdfplot(response(meta2>6.5));
cdfplot(response(meta2<6.5));

%% Graphs of the ActArea

response = trimmedTempDrugData.ActArea;

figure(); hold on
cdfplot(response(meta2>6.5));
cdfplot(response(meta2<6.5));