

T = importDataCCLE();       % reads data
tableT = T;                 % saves a table version for analysis of genes
T = T(:,2:end);             % removes junk
T =  table2array(T);        % converts the table to an array
T = T';                     % formats correctly
[coeff,score,latent] = pca(T);  % performs pca
toPlotT = T*coeff(:,1:3);        % takes three principal components
K=3;                            % makes three clusters
idx = kmeans (toPlotT,K);        % performs k-means clustering
for i = 1:length(idx)     % recombines clusters used to prevent swiss roll
    if (idx(i) == 3)
        idx(i) =2;
    end
    if (idx(i) == 1)
        idx(i) =3;
    end
end

figure(); hold on       % plots clusters using three principal components
for i  = 1:2
    plot3(toPlotT(idx==i+1,1),toPlotT(idx==i+1,2), toPlotT(idx==i+1,3)....
        ,'.', 'MarkerSize',12);
end

xlabel('pca1','FontSize', 10);
ylabel('pca2','FontSize', 10);
zlabel('pca3','FontSize', 10);
title ('Cancer Gene Expression', 'FontWeight', 'bold', 'FontSize',14);


avgExp = mean (T);              % computes the average gene expression
avgExpGroup2 = mean(T(idx==2,:));   %analyzes top up and down regulation
avgExpGroup3 = mean(T(idx==3,:));
names = tableT(:,1);
names = table2array(names);


[values, inds] = sort((avgExp - avgExpGroup2),'descend');
downRegulators2  = names(inds(1:40),1);
upRegulators2 = names(inds((end-39):end),1);

[values, inds] = sort((avgExp - avgExpGroup3),'descend');
downRegulators3  = names(inds(1:40));
upRegulators3 = names(inds((end-39):end));

s