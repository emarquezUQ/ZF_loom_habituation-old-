
%%% this script is to get the rasters of the sound, multisensory and
%%% movements raster plots. 

%%% based on the geting_meansNrasters script, so i can keep the same scale.
%%% see below. 

% figure;imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
% saveas(gcf,'rasterplot_f20_CL4_fasthab.tif');


%%% at the moment I am only doing it for f20


figure;imagesc(ZS_f20(clust_f20_CL7_cleaned.clust_f20_CL7_3_cleaned,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_sound.tif');

figure;imagesc(ZS_f20(idx_multisense_cleaned,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_multisense.tif');


%%% Note: the movement one gives a weird scatter plot cause after cleaning
%%% the idxs in ordered in a different way. In need to sort them the right
%%% way so I can see them by fish. 

figure;imagesc(ZS_f20(idx_rsq_Mov_cleaned,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot


idx_rsq_Mov_cleaned2=zeros(size(idx_rsq_Mov_cleaned));
idx_rsq_Mov_cleaned2=find(ismember(idx_rsq_Mov_cleaned,idx_rsq_Mov));
idx_rsq_Mov_cleaned2=idx_rsq_Mov(idx_rsq_Mov_cleaned2);


%%% not it works     

figure;imagesc(ZS_f20(idx_rsq_Mov_cleaned2,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_mov.tif');

%%%%%%%%%%% WAIT!!! there is an error that i need to
%%% repair. most of the ROIs are sorted properly but it seems that some (like 10%)
%%% not...

figure;scatter(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned2,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned2,2),'.');
hold on;



%%% now it seems to work... i will need to change the raspterplot...
idx_rsq_Mov_cleaned3=zeros(size(idx_rsq_Mov));
idx_rsq_Mov_cleaned3=find(ismember(idx_rsq_Mov,idx_rsq_Mov_cleaned));
idx_rsq_Mov_cleaned3=idx_rsq_Mov(idx_rsq_Mov_cleaned3);

%%% checking the location
figure;scatter(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,2),'.');
hold on;

%%% also checking the numbers
size(intersect(idx_rsq_Mov_cleaned2,idx_rsq_Mov_cleaned)) %%% so this was the wrong one
size(intersect(idx_rsq_Mov_cleaned3,idx_rsq_Mov_cleaned)) %%% and this the correct one





%%

%%% this part is to make the mov raspterplot again as there was a problem boefore and to get the mean of
%%% ohe of the fish. 

load('f20_cleaned_idxs.mat');

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');

load('All_More_BrainReg2.mat');

load('Tail_mov_F20.mat','idx_rsq_Mov');

load('means_F20_CL4n7.mat')


idx_rsq_Mov_cleaned2=zeros(size(idx_rsq_Mov_cleaned));
idx_rsq_Mov_cleaned2=find(ismember(idx_rsq_Mov_cleaned,idx_rsq_Mov));
idx_rsq_Mov_cleaned2=idx_rsq_Mov(idx_rsq_Mov_cleaned2);


%%%%%%%%%%% WAIT!!! there is an error that i need to
%%% repair. most of the ROIs are sorted properly but it seems that some (like 10%)
%%% not...

figure;scatter(ROI_temp2.f20(idx_rsq_Mov_cleaned,1),ROI_temp2.f20(idx_rsq_Mov_cleaned,2),'filled');
hold on;
scatter(ROI_temp2.f20(idx_rsq_Mov_cleaned2,1),ROI_temp2.f20(idx_rsq_Mov_cleaned2,2),'filled');
hold on;


%%% now it seems to work... 
idx_rsq_Mov_cleaned3=zeros(size(idx_rsq_Mov));
idx_rsq_Mov_cleaned3=find(ismember(idx_rsq_Mov,idx_rsq_Mov_cleaned));
idx_rsq_Mov_cleaned3=idx_rsq_Mov(idx_rsq_Mov_cleaned3);


figure;scatter(ROI_temp2.f20(idx_rsq_Mov_cleaned,1),ROI_temp2.f20(idx_rsq_Mov_cleaned,2),'filled');
hold on;
scatter(ROI_temp2.f20(idx_rsq_Mov_cleaned3,1),ROI_temp2.f20(idx_rsq_Mov_cleaned3,2),'filled');
hold on;


figure;imagesc(ZS_f20(idx_rsq_Mov_cleaned3,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_f20_mov_corrected.tif');

%%% now the mean of one of the fish

unique(idx_Fish_f20)

%%% I will try with fish 1 and 44


%figure;plot(mean(ZS_f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),:)));

%%% to get the values to put in prism
mean(ZS_f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==1)),:));

mean(ZS_f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),:));

%%% to see if they are representative
figure;scatter(ROI_temp2.f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==1)),1),ROI_temp2.f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==1)),2),'filled');
hold on;
scatter(ROI_temp2.f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),1),ROI_temp2.f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),2),'filled');
hold on;

