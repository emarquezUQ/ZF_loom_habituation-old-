%%% this is to try to check at the activity in the nMLF per individual fish
%%% to see if I can find motor responses. 

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


ZS_f20=load('f20_CN_r2050_CL5_extra.mat','ZS_CN');
%load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
load('final_F20_step1.mat','idx_Fish_f20');

fish=unique(idx_Fish_f20);



extraBrain=load('All_More_BrainReg4.mat');

extraBrain2=load('All_More_BrainReg3.mat');

%%% it seems that the nMLF could be reg95... although i have an extra one
%%% 145 instead of 144... so I will need to check if I forgot to delete one
%%% and this are in the right location. 


figure;
scatter(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),'.');
hold on;
scatter(Zbrain_Masks{96,3}(:,1),Zbrain_Masks{96,3}(:,2),'.');
hold on;
scatter3(extraBrain.ROI_temp2.f20(extraBrain.PerBrainRegions.f20.reg95.idx,1),extraBrain.ROI_temp2.f20(extraBrain.PerBrainRegions.f20.reg95.idx,2),extraBrain.ROI_temp2.f20(extraBrain.PerBrainRegions.f20.reg95.idx,3),'.');

%%% reg95 is correct
figure;
imagesc(ZS_f20.ZS_CN(extraBrain.PerBrainRegions.f20.reg95.idx,:));

%%
figure;
scatter(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_Mov_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_Mov_cleaned,2),'.');
hold on;
scatter3(extraBrain.ROI_temp2.f20(extraBrain.PerBrainRegions.f20.reg95.idx,1),extraBrain.ROI_temp2.f20(extraBrain.PerBrainRegions.f20.reg95.idx,2),extraBrain.ROI_temp2.f20(extraBrain.PerBrainRegions.f20.reg95.idx,3),'.');

%%% it seems that none of the nmlf ROIs passed the filter of movement
%%% regressors. 
nmlf_mov_idx=intersect(f20_cleaned_idxs.idx_rsq_Mov_cleaned,extraBrain.PerBrainRegions.f20.reg95.idx);


%%% to plot the means of the activity of each fish in the nMLF
counter=1;
figure;
for i=1:length(fish)
    temp_idx=intersect(find(idx_Fish_f20==fish(i)),extraBrain.PerBrainRegions.f20.reg95.idx);
    subplot(3,4,counter);
    plot(mean(ZS_f20.ZS_CN(temp_idx,:)));title(strcat('fish=',num2str(fish(i))));
    
    
    counter=counter+1;
end

%% to look also at hindbrain activity

extraBrain2=load('All_More_BrainReg3.mat');

Rhombs_list=extraBrain2.RegionList(106:111);
regnum=[106:111];

%figure;
fast_hindB_idx_f20=[];
for i=1:length(Rhombs_list)
fast_hindB_idx_f20_temp=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,extraBrain2.PerBrainRegions.f20.(strcat('reg',num2str(regnum(i)))).idx);
fast_hindB_idx_f20=vertcat(fast_hindB_idx_f20,fast_hindB_idx_f20_temp);

% scatter3(extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,1),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,2),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,3),'.');
% hold on;   
    
end
figure;
scatter3(extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20,1),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20,2),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20,3),'.');


%figure;
motor_hindB_idx_f20=[];
for i=1:length(Rhombs_list)
motor_hindB_idx_f20_temp=intersect(f20_cleaned_idxs.idx_rsq_Mov_cleaned,extraBrain2.PerBrainRegions.f20.(strcat('reg',num2str(regnum(i)))).idx);
motor_hindB_idx_f20=vertcat(motor_hindB_idx_f20,motor_hindB_idx_f20_temp);

% scatter3(extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,1),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,2),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,3),'.');
% hold on;   
    
end
figure;
scatter3(extraBrain2.ROI_temp2.f20(motor_hindB_idx_f20,1),extraBrain2.ROI_temp2.f20(motor_hindB_idx_f20,2),extraBrain2.ROI_temp2.f20(motor_hindB_idx_f20,3),'.');


motorR1_hindB_idx_f20=[];

motorR1_hindB_idx_f20_temp=intersect(f20_cleaned_idxs.idx_rsq_Mov_cleaned,extraBrain2.PerBrainRegions.f20.(strcat('reg',num2str(105))).idx);
motorR1_hindB_idx_f20=vertcat(motorR1_hindB_idx_f20,motorR1_hindB_idx_f20_temp);

% scatter3(extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,1),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,2),extraBrain2.ROI_temp2.f20(fast_hindB_idx_f20_temp,3),'.');
% hold on;   
    

figure;
scatter3(extraBrain2.ROI_temp2.f20(motorR1_hindB_idx_f20,1),extraBrain2.ROI_temp2.f20(motorR1_hindB_idx_f20,2),extraBrain2.ROI_temp2.f20(motorR1_hindB_idx_f20,3),'.');


%%
figure;
scatter3(extraBrain2.ROI_temp2.f20(motor_hindB_idx_f20,1),extraBrain2.ROI_temp2.f20(motor_hindB_idx_f20,2),extraBrain2.ROI_temp2.f20(motor_hindB_idx_f20,3),'.');
hold on;
scatter3(extraBrain2.ROI_temp2.f20(motorR1_hindB_idx_f20,1),extraBrain2.ROI_temp2.f20(motorR1_hindB_idx_f20,2),extraBrain2.ROI_temp2.f20(motorR1_hindB_idx_f20,3),'.');

%%
%%% to plot the means of the activity of each fish in the nMLF and the
%%% hindbrain positive for the motor regressors
counter=1;
figure;
for i=1:length(fish)
    subplot(3,4,counter);
    
    temp_idx=intersect(find(idx_Fish_f20==fish(i)),motorR1_hindB_idx_f20);
    plot(mean(ZS_f20.ZS_CN(temp_idx,:)));
    hold on;
    temp_idx=intersect(find(idx_Fish_f20==fish(i)),motor_hindB_idx_f20);
    plot(mean(ZS_f20.ZS_CN(temp_idx,:)));
    hold on;
    temp_idx=intersect(find(idx_Fish_f20==fish(i)),extraBrain.PerBrainRegions.f20.reg95.idx);
    plot(mean(ZS_f20.ZS_CN(temp_idx,:)));title(strcat('fish=',num2str(fish(i))));
       
    
    counter=counter+1;
end


%%%% for some of them it correlates quite nicely!!! i am not sure why they
%%%% didnt got included in the motor cluster... 

%%% but knowing better which regions are responding maybe I can try to
%%% target them in the 60ISI clusters. 