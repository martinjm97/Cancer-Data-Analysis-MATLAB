expression = readtable('CCLE_Expression_Entrez_2012-09-29.gct',...
    'FileType','text',...
    'Delimiter','\t',...
    'HeaderLines',2,...
    'ReadVariableNames',true,...
    'ReadRowNames',true);
%%
names = expression.Description;
names(strcmp('',names)) = expression.Properties.RowNames(strcmp('',names)); 
% replace empty description with probe name
expression.Description = []; % all numeric table
%%
N = 50;
start = zeros(18988,N);
start(1:N,1:N) = eye(N); % this is like 'robust' start for the first N (50) genes
opts = struct('Display','iter'); 
% shows the information for each iteration of metafeatures

%implements the metafeatures algorithm with a reasonable alpha value and 
% random initialization.
[M, W, Gsorted, GsortedInd] = metafeatures...
    (expression,'Alpha',3,'Start',start,'Options',opts); 
