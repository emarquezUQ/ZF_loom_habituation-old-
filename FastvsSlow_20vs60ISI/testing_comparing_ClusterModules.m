
%%%%% this script is to test using a k-means to make 3-d modules of each cluster ROIs, to then use
%%%%% another k-means to compare the 3d modules locations across datasets
%%%%% to see if i can use that as a similarity measure. 


%%% at the moment is not working... but maybe I am doing something wrong. 

load('All_More_BrainReg2.mat');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');



datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);


%%% testing
options = statset('UseParallel',1);
[idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,:),100,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
scatter(Cmap_ROIs(:,1),Cmap_ROIs(:,2));



idxKmeansNCmap=struct;

for i=1:4


idx_clust=find(strcmp(clustersF, clustersF(i)));

temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

[idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),90,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
scatter(Cmap_ROIs(:,1),Cmap_ROIs(:,2));

idxKmeansNCmap.f20.(char(clustersF(i))).idxKmeans_ROIs=idxKmeans_ROIs;
idxKmeansNCmap.f20.(char(clustersF(i))).Cmap_ROIs=Cmap_ROIs;

temp_fieldname=fieldnames(f60_cleaned_idxs.clust_f60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

[idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname)),:),90,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
scatter(Cmap_ROIs(:,1),Cmap_ROIs(:,2));

idxKmeansNCmap.f60.(char(clustersF(i))).idxKmeans_ROIs=idxKmeans_ROIs;
idxKmeansNCmap.f60.(char(clustersF(i))).Cmap_ROIs=Cmap_ROIs;


idx_clust=find(strcmp(clustersS, clustersF(i)));

temp_fieldname=fieldnames(s20_cleaned_idxs.clust_s20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

[idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname)),:),90,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
scatter(Cmap_ROIs(:,1),Cmap_ROIs(:,2));

idxKmeansNCmap.s20.(char(clustersF(i))).idxKmeans_ROIs=idxKmeans_ROIs;
idxKmeansNCmap.s20.(char(clustersF(i))).Cmap_ROIs=Cmap_ROIs;

temp_fieldname=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

[idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname)),:),90,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
scatter(Cmap_ROIs(:,1),Cmap_ROIs(:,2));

idxKmeansNCmap.s60.(char(clustersF(i))).idxKmeans_ROIs=idxKmeans_ROIs;
idxKmeansNCmap.s60.(char(clustersF(i))).Cmap_ROIs=Cmap_ROIs;


end

%%

counter=1;
Cmap_Mat=[];
%Cmap_Mat(1,:) = idxKmeansNCmap.f20.nonhab.Cmap_ROIs(:);


for i=1:4
for j=1:4

    
Cmap_Mat(counter,:) = idxKmeansNCmap.(datasets(i,:)).(char(clustersF(j))).Cmap_ROIs(:);    

counter=counter+1

end


end

%%

%%% to check the location of the modules
counter=1;
figure;
for i=1:4
for j=1:4

subplot(4,4,counter);    
temp_mod=idxKmeansNCmap.(datasets(i,:)).(char(clustersF(j))).Cmap_ROIs;    

scatter(temp_mod(:,1),temp_mod(:,2),'.');

counter=counter+1

end

end


%%
[idxKmeans_ROIs Cmap_ROIs]=kmeans(Cmap_Mat,4,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

mat={};
for i=1:4

    mat{i} = vec2mat(Cmap_ROIs(i,:),3);

end
scatter(mat{1,1}(:,1),mat{1,1}(:,2)); %%%% mmm doesnt work... i was hoping it would make a 'model' of each cluster... 
