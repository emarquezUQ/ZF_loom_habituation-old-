%%% for CL4




load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx', 'S_trim');

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

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
%BrainReg_F20=load('BrainReg_F20.mat');

%PerBrainRegions_f20=BrainReg_F20.PerBrainRegions;


%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_f20_CL4_cleaned=f20_cleaned_idxs.clust_f20_CL4_cleaned;

clust_f20_CL4_cleaned_cell={};
clust=fieldnames(clust_f20_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL4_cleaned_cell.(clustersF{j,1})=clust_f20_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_f20=struct;
Means_CL4_per_fishNbrain_f20=struct;
fish=unique(idx_Fish_f20);


for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f20.(RegionList{brain}).idx);

if length(temp_idx)==1
temp_mean=ZS_f20(temp_idx,:); 
else
temp_mean=mean(ZS_f20(temp_idx,:));
end
Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(tempfish,:)=temp_mean;

%%% just for testing how to do it, first if i want a field per fish and
%%% then if i want them as a matrix in the cluster field.
%test_idx_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}).(strcat('fish',num2str(fish(tempfish))))=temp_idx;
%test_idx_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(tempfish,:)=mean(temp_idx);


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. they
%%% look cool! there are differences as expeted per brain region. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f20{1,clust}))

end
end


%%% now ploting the means.
%%% the results are very interesting too

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f20_perfishNbrain=struct;

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1));
end


Max_resp_f20_perfishNbrain.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Max_resp_f20_perfishNbrain.(RegionList{brain}).(clustersF{clust})(i,:)); 
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
for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
subplot(3,3,counter);
plot(nanmean(Max_resp_f20_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain




%%
%%%for f60

load('final_F60_step1_2.mat','ZS_f60','idx_Fish_f60','ZS_short_F60');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
%BrainReg_F60=load('BrainReg_F60.mat');

%PerBrainRegions_f60=BrainReg_F60.PerBrainRegions;


%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_f60_CL4_cleaned=f60_cleaned_idxs.clust_f60_CL4_cleaned;

clust_f60_CL4_cleaned_cell={};
clust=fieldnames(clust_f60_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f60_CL4_cleaned_cell.(clustersF{j,1})=clust_f60_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_f60=struct;
Means_CL4_per_fishNbrain_f60=struct;
fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...


for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f60_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f60.(RegionList{brain}).idx);

if length(temp_idx)==1
temp_mean=ZS_f60(temp_idx,ZS_short_F60); 
else
temp_mean=mean(ZS_f60(temp_idx,ZS_short_F60));
end
Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(tempfish,:)=temp_mean;

%%% just for testing how to do it, first if i want a field per fish and
%%% then if i want them as a matrix in the cluster field.
%test_idx_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust}).(strcat('fish',num2str(fish(tempfish))))=temp_idx;
%test_idx_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(tempfish,:)=mean(temp_idx);


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. they
%%% look cool! there are differences as expeted per brain region. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f60{1,clust}))

end
end


%%% now ploting the means.
%%% the results are very interesting too

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f60_perfishNbrain=struct;

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1));
end


Max_resp_f60_perfishNbrain.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Max_resp_f60_perfishNbrain.(RegionList{brain}).(clustersF{clust})(i,:)); 
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
for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
subplot(3,3,counter);
plot(nanmean(Max_resp_f60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain



%%

%%% now for s20

load('final_S20_step1.mat','ZS_s20','idx_Fish_s20');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
% BrainReg_S20=load('BrainReg_S20.mat');
% 
% PerBrainRegions_s20=BrainReg_S20.PerBrainRegions;


%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_s20_CL4_cleaned=s20_cleaned_idxs.clust_s20_CL4_cleaned;

clust_s20_CL4_cleaned_cell={};
clust=fieldnames(clust_s20_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s20_CL4_cleaned_cell.(clustersS{j,1})=clust_s20_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_s20=struct;
Means_CL4_per_fishNbrain_s20=struct;
fish=unique(idx_Fish_s20);



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s20_CL4_cleaned_cell.(clustersS{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s20.(RegionList{brain}).idx);

if length(temp_idx)==1
temp_mean=ZS_s20(temp_idx,S_trim); 
else
temp_mean=mean(ZS_s20(temp_idx,S_trim));
end
Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersS{clust})(tempfish,:)=temp_mean;

%%% just for testing how to do it, first if i want a field per fish and
%%% then if i want them as a matrix in the cluster field.
%test_idx_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust}).(strcat('fish',num2str(fish(tempfish))))=temp_idx;
%test_idx_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(tempfish,:)=mean(temp_idx);


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. they
%%% look cool! there are differences as expeted per brain region. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(i,:));  %%% i use clustersF to have it in the same order than f60
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_s20{1,clust}))

end
end


%%% now ploting the means.
%%% the results are very interesting too

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s20_perfishNbrain=struct;

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersS{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(k)+1));
end


Max_resp_s20_perfishNbrain.(RegionList{brain}).(clustersS{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
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
   plot(Max_resp_s20_perfishNbrain.(RegionList{brain}).(clustersF{clust})(i,:)); %%% i use clustersF to have it in the same order than f60
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
plot(nanmean(Max_resp_s20_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain Means_CL4_per_fishNbrain_s20 Max_resp_s20_perfishNbrain

%%


%%% now for s60

load('final_S60_step1.mat','ZS_s60','idx_Fish_s60','ZS_short_S60');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');
% BrainReg_S60=load('BrainReg_S60.mat');
% 
% PerBrainRegions_s60=BrainReg_S60.PerBrainRegions;



%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_s60_CL4_cleaned=s60_cleaned_idxs.clust_s60_CL4_cleaned;

clust_s60_CL4_cleaned_cell={};
clust=fieldnames(clust_s60_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s60_CL4_cleaned_cell.(clustersS{j,1})=clust_s60_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_s60=struct;
Means_CL4_per_fishNbrain_s60=struct;
fish=unique(idx_Fish_s60);



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s60_CL4_cleaned_cell.(clustersS{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s60.(RegionList{brain}).idx);

if length(temp_idx)==1
temp_mean=ZS_s60(temp_idx,ZS_short_S60(S_trim)); 
else
temp_mean=mean(ZS_s60(temp_idx,ZS_short_S60(S_trim)));
end
Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersS{clust})(tempfish,:)=temp_mean;

%%% just for testing how to do it, first if i want a field per fish and
%%% then if i want them as a matrix in the cluster field.
%test_idx_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust}).(strcat('fish',num2str(fish(tempfish))))=temp_idx;
%test_idx_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(tempfish,:)=mean(temp_idx);


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. they
%%% look cool! there are differences as expeted per brain region. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(i,:));  %%% i use clustersF to have it in the same order than f60
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_s60{1,clust}))

end
end


%%% now ploting the means.
%%% the results are very interesting too

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
subplot(3,3,counter);
plot(nanmean(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s60_perfishNbrain=struct;

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersS{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersS{clust})(i,Loomf20_onset_idx(k)+1));
end


Max_resp_s60_perfishNbrain.(RegionList{brain}).(clustersS{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
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
   plot(Max_resp_s60_perfishNbrain.(RegionList{brain}).(clustersF{clust})(i,:)); %%% i use clustersF to have it in the same order than f60
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
plot(nanmean(Max_resp_s60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain Means_CL4_per_fishNbrain_s20 Max_resp_s20_perfishNbrain Means_CL4_per_fishNbrain_s60 Max_resp_s60_perfishNbrain


%%



save('Max_response_perBrain_tryingthings_all_CL4.mat');



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
plot(nanmean(Max_resp_f60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s20_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


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
plot(nanmean(Max_resp_f20_perfishNbrain.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_f60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s20_perfishNbrain.(RegionList{brain}).(clustersF{clust})));
% hold on;
% plot(nanmean(Max_resp_s60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end
end


%%% now seeing it by figures for each cluster and then by brain region. 

mkdir clustersPerBrainRegion_CL4
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


subplot(3,3,counter);
for fish=1:size(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}),1)
plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(fish,:)); title(strcat('f20_',RegionList{brain}));
hold on;
end
% 


% subplot(3,3,counter);
% for fish=1:size(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust}),1)
% plot(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(fish,:)); title(strcat('f60_',RegionList{brain}));
% hold on;
% end
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
plot(nanmean(Max_resp_f20_perfishNbrain.(RegionList{brain}).(clustersF{clust}))); title(RegionList{brain});
hold on;
plot(nanmean(Max_resp_f60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s20_perfishNbrain.(RegionList{brain}).(clustersF{clust})));
hold on;
plot(nanmean(Max_resp_s60_perfishNbrain.(RegionList{brain}).(clustersF{clust})));


counter=counter+1;


end

%saveas(gcf,strcat(savedir,'\MaxPerGroup_',clustersF{clust}),'png');
end
close all
