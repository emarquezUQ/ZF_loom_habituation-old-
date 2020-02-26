%%%% this script is to get calculate correlations between the tectum Max
%%%% responses and the behavioural responeses from the free swiming
%%%% animals. this time is to do it per fish so I can calculate confidence
%%%% intervals and be able to conclude that the weakly habituating one is
%%%% the cluster that most looks like the behavioural responses. 

%%% I will also try to do this per brain region. 

load('All_means_maxResp_CL4_normalized.mat');

%%%% reordering the RegionList

RegionList2=RegionList([1 2 3 9 4 5 8 7 6]);
RegionList=RegionList2;
%%%% get the data from the behaviour (I took it from prism and made a new
%%%% variable)

Behaviour=[];

%%%% then I put in the variable the values from the Prism file
%%%% (FnS_20vs60ISI_20191208), per columns, f20,f60,s20,s60


Correlations_BrainVSbehav=struct;

for brain=1:length(RegionList)

counter=1;


for clust=1:3

if ismember(clustersF{clust},fieldnames(Max_resp_f20_perfishNbrain2.(RegionList{brain})))
 for f=1:size(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.f20.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
 end
else
    0 %%%% the zeros are just to flag that there are some cases where I wont have data.   
end 

if ismember(clustersF{clust},fieldnames(Max_resp_f60_perfishNbrain2.(RegionList{brain})))
for f=1:size(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.f60.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
end
else
    0
end

if ismember(clustersF{clust},fieldnames(Max_resp_s20_perfishNbrain2.(RegionList{brain})))
for f=1:size(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.s20.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
end
else
    0
end

if ismember(clustersF{clust},fieldnames(Max_resp_s60_perfishNbrain2.(RegionList{brain})))
for f=1:size(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.s60.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
end
else
    0
end

counter=counter+1;

end

end
clear temp

%%%% visualizing the results
%%% I am asking that there are at least 3 fish
Corr_BrainVSbehav_mat=struct;
Corr_mat_mean=[];
for data=1:length(datasets)
    
    for brain=1:length(RegionList)
        counter=1;
        for clust=[2 3 1] %%% to put the clusters in order (fasthab, slopehab, nonhab).
            if ismember(clustersF{clust},fieldnames(Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain})))
            if 3<sum(~isnan(Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain}).(clustersF{clust})(:,1)));
            temp=nanmean(Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain}).(clustersF{clust})(:,1));
            else
            temp=NaN;
            end
            else
            temp=NaN;
            end
            Corr_BrainVSbehav_mat.(datasets(data,:))(counter,brain)=temp;
            counter=counter+1;
        end
    end
    
    Corr_mat_mean=cat(3,Corr_mat_mean,Corr_BrainVSbehav_mat.(datasets(data,:)));
    Corr_mat_mean=mean(Corr_mat_mean,3); %%%% note, as I am doing a normal mean I am leaving out the combinations that have NaNs. so this means that for it to be included there most be at least 3 fish in the 4 datasets. 
    
    figure;imagesc(Corr_BrainVSbehav_mat.(datasets(data,:))),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));
end
clear temp
figure;imagesc(Corr_mat_mean),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));

%%%% now with inferno colormap. need to add to path 
figure;imagesc(Corr_mat_mean),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));colorbar;colormap(inferno);
saveas(gcf,'Max_resp_Brain_corr_freeswim.svg');

%%%%% what if I treat NaN as 0s? but I dont think we can use this

%%%% visualizing the results
% Corr_BrainVSbehav_mat=struct;
% Corr_mat_mean=[];
% for data=1:length(datasets)
%     
%     for brain=1:length(RegionList)
%         counter=1;
%         for clust=[2 3 1] %%% to put the clusters in order (fasthab, slopehab, nonhab).
%             if ismember(clustersF{clust},fieldnames(Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain})))
%             temp1=Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain}).(clustersF{clust})(:,1);
%             temp1(find(isnan(temp1)))=0;    
%             temp=mean(temp1);          
%             else
%             temp=NaN;
%             end
%             Corr_BrainVSbehav_mat.(datasets(data,:))(counter,brain)=temp;
%             counter=counter+1;
%             clear temp temp1
%         end
%     end
%     Corr_mat_mean=cat(3,Corr_mat_mean,Corr_BrainVSbehav_mat.(datasets(data,:)));
%     Corr_mat_mean=mean(Corr_mat_mean,3);
%     figure;imagesc(Corr_BrainVSbehav_mat.(datasets(data,:))),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));
%     
% end
% figure;imagesc(Corr_mat_mean),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));
% clear temp temp1