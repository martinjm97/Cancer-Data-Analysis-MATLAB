
T = importData; % imports data from various files and returns a table
tableT = T;     % stores the table for later use
T =  table2array(T);        % converts the table to an array
T = knnimpute(T);           % fills in missing data
T = T';              % puts data in form (num of samples) x (num of genes)
[coeff,score,latent] = pca(T);  % performs pca
toPlot = T*coeff(:,1:3);        % takes three principal components
K=3;                            % makes three clusters 
idx = kmeans (toPlot,K);        % performs k-means clustering

avgExp = mean (T);              % computes the average gene expression

% computes the average gene expression for each group
avgExpGroup1 = mean(T(idx==1,:));
avgExpGroup2 = mean(T(idx==2,:));
avgExpGroup3 = mean(T(idx==3,:));

% plots data using two principal components and three clusters
figure(); hold on 
for i  = 1:3
    plot(toPlot(idx==i,1),toPlot(idx==i,2),'.'...
        ,'MarkerSize',12);
end

xlabel('pca1','FontSize', 10);
ylabel('pca2','FontSize', 10);
zlabel('pca3','FontSize', 10);
title ('Breast Cancer Gene Expression', 'FontWeight', 'bold', 'FontSize',14);



isRegulator1 = (abs(avgExp - avgExpGroup1) > 1);

consecNums = [1: length(isRegulator1)];

% lists top fourty up and down regulated genes
regulators1 = tableT.Properties.RowNames(consecNums (isRegulator1));

isRegulator2 = (abs(avgExp - avgExpGroup2) > 3);
regulators2 = tableT.Properties.RowNames(consecNums (isRegulator2));

isRegulator3 = (abs(avgExp - avgExpGroup3) > 3);
regulators3 = tableT.Properties.RowNames(consecNums (isRegulator3));