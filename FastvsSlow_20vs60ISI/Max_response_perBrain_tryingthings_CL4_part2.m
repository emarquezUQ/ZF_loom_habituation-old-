%%% this script is a follow up of the Max_response_perBrain_tryingthings_CL4.m and is just to
%%% add a timepoint 0 for my max response graphs. I am taking the onset + 1
%%% of the first loom. 

%%% for CL4




load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx', 'S_trim');

load('Max_response_perBrain_tryingthings_all_CL4.mat');



datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);

load('Zbrain_Masks.mat');
load('All_More_BrainReg.mat','PerBrainRegions');


%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain'};

RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};

%%
%%%for f20

load('final_F20_step1.mat','idx_Fish_f20');


fish=unique(idx_Fish_f20);


%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f20_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersF)

for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1));

for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint zero
end


Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);%%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
for i=1:length(fish)
   plot(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
end


counter=counter+1;
end
end


%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
plot(nanmean(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


%clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain




%%
%%%for f60

load('final_F60_step1_2.mat','idx_Fish_f60','ZS_short_F60');

fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...




%%% now ploting the means.
%%% the results are very interesting too

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f60_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersF)

for  i=1:length(fish) 
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1));

    
for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1));
end


Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);
end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
for i=1:length(fish)
   plot(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
end


counter=counter+1;
end
end


%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
plot(nanmean(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


%clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain



%%

%%% now for s20

load('final_S20_step1.mat','idx_Fish_s20');

fish=unique(idx_Fish_s20);



%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s20_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersS)

for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(1)+1));

for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersS{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(k)+1));
end


Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersS{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);
end
end
end




%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
plot(nanmean(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

%clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain Means_CL4_per_fishNbrain_s20 Max_resp_s20_perfishNbrain

%%


%%% now for s60

load('final_S60_step1.mat','idx_Fish_s60','ZS_short_S60');

fish=unique(idx_Fish_s60);




%%% now ploting the means.
%%% the results are very interesting too

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersS)
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s60_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersS)

for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(1)+1));

    
for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersS{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(k)+1));
end


Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersS{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);
end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
for i=1:length(fish)
   plot(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)); %%% i use clustersF to have it in the same order than f60
    hold on;
end


counter=counter+1;
end
end


%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
plot(nanmean(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain2 Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain2 Means_CL4_per_fishNbrain_s20 Max_resp_s20_perfishNbrain2 Means_CL4_per_fishNbrain_s60 Max_resp_s60_perfishNbrain2


%%



save('Max_response_perBrain_tryingthings_all_CL4_part2.mat');



%%

%%% now to compare f60 vs s20


for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end





for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(3,3,counter);
plot(nanmean(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end



%%

%%% now to compare all
 %%% means of each fish. 
for brain=1:length(RegionList)

for clust=1:length(clustersF)
figure('Position',[100 0 900 900]);suptitle(strcat(RegionList{brain},'_',clustersF{clust}));
counter=1;

subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(fish,:)); title('f20');
hold on;
end

counter=counter+1;

subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(fish,:)); title('f60');
hold on;
end

counter=counter+1;


subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(fish,:)); title('s20');
hold on;
end

counter=counter+1;


subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(fish,:)); title('s60');
hold on;
end



end
end


%%% means per cluster
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(2,2,counter);
plot(nanmean(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end





for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(clustersF)
subplot(2,2,counter);
plot(nanmean(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));
% hold on;
% plot(nanmean(Max_resp_s60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


%%% now seeing it by figures for each cluster and then by brain region. 

% mkdir clustersPerBrainRegion_CL4
savedir=('Y:\Emmanuel_MeDiCi\FvsS_20vs60_CNMF\matlab\clustersPerBrainRegion_CL4');
%%%

%%% means of each fish. 
for clust=1:length(clustersF)

for  brain=1:length(RegionList)
figure('Position',[100 0 900 900]);suptitle(strcat(clustersF{clust},'_',RegionList{brain}));
counter=1;

subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(fish,:)); title('f20');
hold on;
end

counter=counter+1;

subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(fish,:)); title('f60');
hold on;
end

counter=counter+1;


subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(fish,:)); title('s20');
hold on;
end

counter=counter+1;


subplot(2,2,counter);
for fish=1:size(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(fish,:)); title('s60');
hold on;
end



end
end



%%% means of each fish in per data set.

for clust=1:length(clustersF)
figure('Position',[100 0 900 900]);suptitle(clustersF{clust});
counter=1;

for  brain=1:length(RegionList)


% subplot(3,3,counter);
% for fish=1:size(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}),1)
% plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(fish,:)); title(strcat('f20_',RegionList{brain}));
% hold on;
% end
% 


subplot(3,3,counter);
for fish=1:size(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(fish,:)); title(strcat('f60_',RegionList{brain}));
hold on;
end
% 


% 
% subplot(3,3,counter);
% for fish=1:size(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust}),1)
% plot(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(fish,:)); title(strcat('s20_',RegionList{brain}));
% hold on;
% end
% 


% subplot(3,3,counter);
% for fish=1:size(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust}),1)
% plot(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(fish,:)); title(strcat('s60_',RegionList{brain}));
% hold on;
% end





counter=counter+1;

end

%saveas(gcf,strcat(savedir,'\MeanPerFish_f20_',clustersF{clust}),'png');
end


%%% means per group
for clust=1:length(clustersF)

figure('Position',[100 0 900 900]);suptitle(clustersF{clust});
counter=1;
for  brain=1:length(RegionList)
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}))); title(RegionList{brain});
hold on;
plot(nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end

%saveas(gcf,strcat(savedir,'\meansPerGroup_',clustersF{clust}),'png');
end

close all

%%% max response

for clust=1:length(clustersF) 

figure('Position',[100 0 900 900]); suptitle(clustersF{clust});
counter=1;
for brain=1:length(RegionList)
subplot(3,3,counter); 
plot(nanmean(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust}))); title(RegionList{brain});
hold on;
plot(nanmean(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end

%saveas(gcf,strcat(savedir,'\MaxPerGroup_',clustersF{clust}),'png');
end
close all

