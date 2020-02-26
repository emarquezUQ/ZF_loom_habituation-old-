%%%%  this is to look at how the clusters look with in each dataset but
%%%%  across individual fish. is a continuation of the
%%%%  testing_ROIsPlots_4datasets.m

load('All_More_BrainReg2.mat');
load('specialMasks.mat');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');


datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);




%%% in f20, i then noticed that the fasthab cluster neurons in the hind brain are
%%% only located in some of the fish (4 out of 11). it seems to be the ones that had an
%%% escape response to the first loom! except for fish 40... which is weird.
%%% I checked its movie and indeed it doesnt move at the first loom. 

%%%I will try to confirm that with other datasets... mmm it seems that in
%%%s20 the ROIs in the hindbrain are more common. and in f60 is less common
%%%but they dont seem to correlate as well with the movement to the fist
%%%loom as in f20

ZS_f20=load('f20_CN_r2050_CL5_extra.mat','ZS_CN');
%load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
load('final_F20_step1.mat','idx_Fish_f20');


fish=unique(idx_Fish_f20);

for i=1:4

figure;
for j=1:length(fish)

   
idx_clust=find(strcmp(clustersF, clustersF(i)));
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

idx_fish_temp=find(idx_Fish_f20==fish(j));

idx_clust_temp=f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

scatter(ROI_temp2.f20(idx_temp,1),ROI_temp2.f20(idx_temp,2),20,'filled');
%scatter3(ROI_temp2.f20(idx_temp,1),ROI_temp2.f20(idx_temp,2),ROI_temp2.f20(idx_temp,3),20,'filled','MarkerFaceAlpha',.2);
hold on;


end

end

%%


for i=1:4

figure;
counter=1;
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersF, clustersF(i)));
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

idx_fish_temp=find(idx_Fish_f20==fish(j));

idx_clust_temp=f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

scatter(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),20,'filled');
hold on;
scatter(ROI_temp2.f20(idx_temp,1),ROI_temp2.f20(idx_temp,2),20,'filled');
%scatter3(ROI_temp2.f20(idx_temp,1),ROI_temp2.f20(idx_temp,2),ROI_temp2.f20(idx_temp,3),20,'filled','MarkerFaceAlpha',.2);
%hold on;

counter=counter+1;
end

end

%%%% there is something interesting... only 4 fish have actually ROIs in the hind
%%%% brain... I will try to see the traces and the movements
load('Tail_mov_F20.mat','Fish_with_behav_f20','Movements');

figure;
counter=1;
i=2; %%% for fasthab
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersF, clustersF(i)));
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

idx_fish_temp=find(idx_Fish_f20==fish(j));

idx_clust_temp=f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

idx_temp=intersect(idx_temp,PerBrainRegions2.hindB_mask.idx);


plot(mean(ZS_f20.ZS_CN(idx_temp,:)));
hold on;
plot(squeeze(squeeze(Movements(find(Fish_with_behav_f20==fish(j)),:))));
hold on;
plot(mean(ZS_f20.ZS_CN(intersect(idx_fish_temp,f20_cleaned_idxs.idx_rsq_Mov_cleaned),:)));

counter=counter+1;
end

%% checking now with s20

load('Tail_mov_S20.mat','Fish_with_behav_s20','Movements'); %%% there is an error when loading... doesnt make any sense


%ZS_s20=load('final_S20_step1.mat','ZS_s20');
load('final_S20_step1.mat','idx_Fish_s20');
S_trim=[1:448 453:901 906:1352];

fish=unique(idx_Fish_s20);


for i=1:4

figure;
counter=1;
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersS, clustersF(i)));
temp_fieldname=fieldnames(s20_cleaned_idxs.clust_s20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);


idx_fish_temp=find(idx_Fish_s20==fish(j));

idx_clust_temp=s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

scatter(ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,1),ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,2),20,'filled');
hold on;
scatter(ROI_temp2.s20(idx_temp,1),ROI_temp2.s20(idx_temp,2),20,'filled');
%scatter3(ROI_temp2.s20(idx_temp,1),ROI_temp2.s20(idx_temp,2),ROI_temp2.s20(idx_temp,3),20,'filled','MarkerFaceAlpha',.2);
%hold on;

counter=counter+1;
end

end

%%%% to plot the movements
figure;
counter=1;
i=2; %%% for fasthab
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersS, clustersF(i)));
temp_fieldname=fieldnames(s20_cleaned_idxs.clust_s20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);


idx_fish_temp=find(idx_Fish_s20==fish(j));

idx_clust_temp=s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

%idx_temp=intersect(idx_temp,PerBrainRegions.hindB_mask.idx);


%plot(mean(ZS_f20.ZS_CN(idx_temp,:)));
%hold on;
plot(squeeze(squeeze(Movements(find(Fish_with_behav_s20==fish(j)),:))));
hold on;
%plot(mean(ZS_f20.ZS_CN(intersect(idx_fish_temp,f20_cleaned_idxs.idx_rsq_Mov_cleaned),:)));

counter=counter+1;
end

%% with f60


load('Tail_mov_F60.mat','Fish_with_behav_f60','Movements');

%ZS_f60=load('final_F60_step1_2.mat','ZS_f60');
load('final_F60_step1_2.mat','ZS_short_F60','idx_Fish_f60');


fish=unique(idx_Fish_f60);

for i=1:4

figure;
counter=1;
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersS, clustersF(i)));
temp_fieldname=fieldnames(f60_cleaned_idxs.clust_f60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);


idx_fish_temp=find(idx_Fish_f60==fish(j));

idx_clust_temp=f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

scatter(ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,1),ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,2),20,'filled');
hold on;
scatter(ROI_temp2.f60(idx_temp,1),ROI_temp2.f60(idx_temp,2),20,'filled');
%scatter3(ROI_temp2.f60(idx_temp,1),ROI_temp2.f60(idx_temp,2),ROI_temp2.f60(idx_temp,3),20,'filled','MarkerFaceAlpha',.2);
%hold on;

counter=counter+1;
end

end

%%%% to plot the movements
figure;
counter=1;
i=2; %%% for fasthab
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersF, clustersF(i)));
temp_fieldname=fieldnames(f60_cleaned_idxs.clust_f60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);


idx_fish_temp=find(idx_Fish_f60==fish(j));

idx_clust_temp=f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

%idx_temp=intersect(idx_temp,PerBrainRegions2.hindB_mask.idx);


plot(mean_CL4_f60.fasthab);
hold on;
plot(squeeze(squeeze(Movements(find(Fish_with_behav_f60==fish(j))',:))));
hold on;
%plot(mean(ZS_f20.ZS_CN(intersect(idx_fish_temp,f20_cleaned_idxs.idx_rsq_Mov_cleaned),:)));

counter=counter+1;
end


%% with s60


load('Tail_mov_S60.mat','Fish_with_behav_s60','Movements');

%ZS_s60=load('final_S60_step1.mat','ZS_s60');
load('final_s60_step1.mat','ZS_short_S60','idx_Fish_s60');


fish=unique(idx_Fish_s60);

for i=1:4

figure;
counter=1;
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersS, clustersF(i)));
temp_fieldname=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);


idx_fish_temp=find(idx_Fish_s60==fish(j));

idx_clust_temp=s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

scatter(ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,1),ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,2),20,'filled');
hold on;
scatter(ROI_temp2.s60(idx_temp,1),ROI_temp2.s60(idx_temp,2),20,'filled');
%scatter3(ROI_temp2.s60(idx_temp,1),ROI_temp2.s60(idx_temp,2),ROI_temp2.s60(idx_temp,3),20,'filled','MarkerFaceAlpha',.2);
%hold on;

counter=counter+1;
end

end

%%%% to plot the movements
figure;
counter=1;
i=2; %%% for fasthab
for j=1:length(fish)
subplot(3,4,counter)
   
idx_clust=find(strcmp(clustersS, clustersF(i)));
temp_fieldname=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);


idx_fish_temp=find(idx_Fish_s60==fish(j));

idx_clust_temp=s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname));

idx_temp=intersect(idx_fish_temp,idx_clust_temp);

%idx_temp=intersect(idx_temp,PerBrainRegions2.hindB_mask.idx);


plot(mean_CL4_s60.fasthab);
hold on;
plot(squeeze(squeeze(Movements(find(Fish_with_behav_s60==fish(j)),:))));
%hold on;
%plot(mean(ZS_f20.ZS_CN(intersect(idx_fish_temp,f20_cleaned_idxs.idx_rsq_Mov_cleaned),:)));

counter=counter+1;
end





%%

%%% i took this bit from testing_ROIsPlots_4datasets.m to use as a template
%%% above. 

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

