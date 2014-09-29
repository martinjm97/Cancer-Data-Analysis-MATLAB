function T = importData ()

% importData() organizes the data contained within the same file.
% returns T which is a table of all of the data.


 dataNames = dir ('*.txt');         % makes array of all of the data files

T = array2table( zeros( 17814,length(dataNames) ) );    % holds data

for iter = 1: length(dataNames)     % iterates through all files
    name = dataNames(iter).name;
    
    temp = readtable(name,...
        'ReadVariableNames', false,...allows for manual naming of variables
        'ReadRowNames', true,...     names rows
        'HeaderLines', 2,...         removes top two junk lines
        'Delimiter','\t',...         tells how to read data
        'Format','%s%f',...          makes sure variables are correct type
        'TreatAsEmpty','null');    % fixes string error in table
    
    T(:,iter) =  temp;             % assigns a column
    T .Properties.VariableNames(iter) = {name(1:23)};   
    % gives a name as part of the file
    
end
T .Properties.RowNames = temp.Properties.RowNames; %names rows correctly
 
end


