%%%%% this script is to get the means of the clusters per brain region to
%%%%% feed the CCM that I will test in R. 


%%% for CL4

load('Max_response_perBrain_all_CL4_corrected.mat')



%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};


%%
%%% for all datasets


%%%% getting the means of each cluster per brain region to feed the CCM
%%%% that I will test in R.

Means_clust_reg_CL4=struct;

for data=1:length(datasets)    
for clust=1:length(clustersF)
figure('Position',[100 0 900 900]);suptitle(clustersF{clust});
counter=1;
for  brain=1:length(RegionList)
subplot(3,3,counter);
if data==1
if ismember(clustersF{clust},fieldnames(Max_resp_f20_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}))); title(RegionList{brain});
Means_clust_reg_CL4.(datasets(data,:)).(clustersF{clust})(:,brain)=nanmean(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust}));
else
end
elseif data==2

if ismember(clustersF{clust},fieldnames(Max_resp_f60_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})));
Means_clust_reg_CL4.(datasets(data,:)).(clustersF{clust})(:,brain)=nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust}));
else
end
elseif data==3
if ismember(clustersF{clust},fieldnames(Max_resp_s20_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})));
Means_clust_reg_CL4.(datasets(data,:)).(clustersF{clust})(:,brain)=nanmean(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust}));
else
end
else
if ismember(clustersF{clust},fieldnames(Max_resp_s60_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})));
Means_clust_reg_CL4.(datasets(data,:)).(clustersF{clust})(:,brain)=nanmean(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust}));
else
end
end
counter=counter+1;


end


end
end
