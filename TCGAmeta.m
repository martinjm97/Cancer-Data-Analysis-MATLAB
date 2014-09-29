%% Import Data

T = importData; % imports data from various files and returns a table


%% Store Data

tableT = T;     % stores the table for later use
geneNames = tableT.Properties.RowNames;
T =  table2array(T);        % converts the table to an array
T = knnimpute(T);           % fills in missing data

%% Performs Metafeatures Algorithm

N = 50;
start = zeros(17814,N);
start(1:N,1:N) = eye(N); % this is like 'robust' start for the first N (50) genes
opts = struct('Display','iter'); % shows the information for each iteration of metafeatures

%implements the metafeatures algorithm with a reasonable alpha value and 
% random initialization.
[M, W, Gsorted, GsortedInd] = metafeatures...
    (T,'Alpha',3,'Start',start,'Options',opts); 

meta = W'*T;   %makes metafeatures

%% Replicates Demo
opts = struct('Display','iter');

erbb         = find(strcmp('ERBB2',geneNames));
estrogen     = find(strcmp('ESR1',geneNames));
progestrone  = find(strcmp('PGR',geneNames));

startValues = eye(size(T,1),3);
startValues(erbb,1)        = 1;
startValues(estrogen,2)    = 1;
startValues(progestrone,3) = 1;

[metaf, weights, genes_sorted] = metafeatures(T, geneNames,'Alpha',5, 'start',startValues,'options',opts);
mets = weights'*T;  

%% Plots from Modified Demo
% compares HER2, ESR1, and ribosomal protien metagenes

figure(); hold on
plot3(metaf(1,:),metaf(2,:),metaf(3,:),'.','MarkerSize',12)
xlabel('ERBB2 metagene')
ylabel('Estrogen metagene')
zlabel('Progesterone')
title ('Breast Cancer Metafeatures', 'FontWeight', 'bold', 'FontSize',14);


% outputs the top 10 genes in each metagene
genes_sorted(1:10,1)
genes_sorted(1:10,2)
genes_sorted(1:10,3)
%% Modified Demo
% Uses only two initialization to find a third cluster

opts = struct('Display','iter');

erbb         = find(strcmp('ERBB2',geneNames));
estrogen     = find(strcmp('ESR1',geneNames));

startValues = eye(size(T,1),3);
startValues(erbb,1)        = 1;
startValues(estrogen,2)    = 1;

[metaf, weights, genes_sorted] = metafeatures(T, geneNames,'Alpha',6, 'start',startValues,'options',opts);
met = weights'*T;  
%% Plots from Modified Demo
% compares HER2, ESR1, and ribosomal protien metagenes

figure(); hold on
plot3(met(1,:),met(2,:),met(3,:),'.','MarkerSize',12)
xlabel('ERBB2 metagene')
ylabel('Estrogen metagene')
zlabel('Ribosomal Protein metagene')
title ('Breast Cancer Metafeatures', 'FontWeight', 'bold', 'FontSize',14);


% outputs the top 10 genes in each metagene
genes_sorted(1:10,1)
genes_sorted(1:10,2)
genes_sorted(1:10,3)

%% Plots from Modified Demo
% Graphs HER2 metagene against the Ribosomal Protien metagene
% ESR1 was ignored because of lack of correlation

figure(); hold on
plot(met(1,:),met(3,:),'.','MarkerSize',12)
xlabel('ERBB2 metagene')
ylabel('Ribosomal Protein metagene')
title ('Breast Cancer Metafeatures', 'FontWeight', 'bold', 'FontSize',14);

% outputs the top 10 genes in each metagene
genes_sorted(1:10,1)
genes_sorted(1:10,2)
genes_sorted(1:10,3)