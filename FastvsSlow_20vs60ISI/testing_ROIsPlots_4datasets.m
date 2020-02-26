%%% this script is to visualize the distribution of ROIs of the different
%%% clusters across the 4 datasets. 


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


figure;
counter=1;
for i=1:4
subplot(2,2,counter);

idx_clust=find(strcmp(clustersF, clustersF(i)));

temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

temp_fieldname=fieldnames(f60_cleaned_idxs.clust_f60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

idx_clust=find(strcmp(clustersS, clustersF(i)));

temp_fieldname=fieldnames(s20_cleaned_idxs.clust_s20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

temp_fieldname=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

counter=counter+1;

end

%%

figure;
counter=1;
for i=1:4


idx_clust=find(strcmp(clustersF, clustersF(i)));

temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,4,counter);
scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);

counter=counter+1;

temp_fieldname=fieldnames(f60_cleaned_idxs.clust_f60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,4,counter);
scatter(ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
counter=counter+1;

idx_clust=find(strcmp(clustersS, clustersF(i)));

temp_fieldname=fieldnames(s20_cleaned_idxs.clust_s20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,4,counter);
scatter(ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
counter=counter+1;

temp_fieldname=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,4,counter);
scatter(ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
counter=counter+1;

%counter=counter+1;

end

%% 

%%%% if I include sound and fasthab subtypes

clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);



figure;
counter=1;
for i=1:7
subplot(2,4,counter);

idx_clust=find(strcmp(clustersF, clustersF(i)));

temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL7_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

temp_fieldname=fieldnames(f60_cleaned_idxs.clust_f60_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL7_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

idx_clust=find(strcmp(clustersS, clustersF(i)));

temp_fieldname=fieldnames(s20_cleaned_idxs.clust_s20_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL7_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

temp_fieldname=fieldnames(s60_cleaned_idxs.clust_s60_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
scatter(ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL7_cleaned.(char(temp_fieldname)),2),20,'filled','MarkerFaceAlpha',.2);
hold on;

counter=counter+1;

end

%%

%%%% and motor

figure;
counter=1;
counter2=8;
for i=1:7


idx_clust=find(strcmp(clustersF, clustersF(i)));

temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,8,counter);
scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL7_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;

temp_fieldname=fieldnames(f60_cleaned_idxs.clust_f60_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,8,counter+counter2);
scatter(ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL7_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;

idx_clust=find(strcmp(clustersS, clustersF(i)));

temp_fieldname=fieldnames(s20_cleaned_idxs.clust_s20_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,8,counter+counter2*2);
scatter(ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL7_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;

temp_fieldname=fieldnames(s60_cleaned_idxs.clust_s60_CL7_cleaned);temp_fieldname=temp_fieldname(idx_clust);
subplot(4,8,counter+counter2*3);
scatter(ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL7_cleaned.(char(temp_fieldname)),1),ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL7_cleaned.(char(temp_fieldname)),2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;


counter=counter+1;

end

counter=8;
counter2=8;
for i=1:4

subplot(4,8,counter);
scatter(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_Mov_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_Mov_cleaned,2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;

subplot(4,8,counter+counter2);
scatter(ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_Mov_cleaned,1),ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_Mov_cleaned,2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;

subplot(4,8,counter+counter2*2);
scatter(ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_Mov_cleaned,1),ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_Mov_cleaned,2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;

subplot(4,8,counter+counter2*3);
scatter(ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_Mov_cleaned,1),ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_Mov_cleaned,2),4,'filled','MarkerFaceAlpha',.2); xlim([400 1300]); ylim([0 600]);
%counter=counter+1;


%counter=counter+1;

end
