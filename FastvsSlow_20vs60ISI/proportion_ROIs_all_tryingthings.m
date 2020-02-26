
%%% this script is to look at the number/proportion of ROIs across fish and
%%% datasets and to be able to look also at in across brain regions. 

%%%% first you need to load what you need

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab


datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

%load('Zbrain_Masks.mat');

RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain'};

%%
%%% for f20


load('final_F20_step1.mat','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
BrainReg_F20=load('BrainReg_F20.mat');

PerBrainRegions_f20=BrainReg_F20.PerBrainRegions;

idx_rsq_test_f20short_cleaned=f20_cleaned_idxs.idx_rsq_test_f20short_cleaned;
idx_f20_multisense_cleaned=f20_cleaned_idxs.idx_multisense_cleaned;
idx_rsq_Mov_cleaned=f20_cleaned_idxs.idx_rsq_Mov_cleaned;
clust_f20_CL7_cleaned=f20_cleaned_idxs.clust_f20_CL7_cleaned;

clust_f20_CL7_cleaned_cell={};
clust=fieldnames(clust_f20_CL7_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL7_cleaned_cell.(clustersF{j,1})=clust_f20_CL7_cleaned.(clust{j});   
end    



%%% to look at the proportion per brain region to see if it works
%%% NOTE: i am basing the proportion on the response rsq filtered ROIs, not in
%%% the total number of ROIs in that brian region. it could be intersting
%%% to look at both.

counter=1;
figure;
for j=1:size(clustersF,1)
    
    idx_temp=clust_f20_CL7_cleaned_cell.(clustersF{j,1});
    
    for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_f20.(regionName).idx,idx_temp))/sum(ismember(PerBrainRegions_f20.(regionName).idx,idx_rsq_test_f20short_cleaned));
end
subplot(3,3,counter);bar(ratioperRegion);set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{j,1})); ylim([0 1])
    
   %hold on; 
   counter=counter+1;
end

%%% now i need to do it an average per fish

Prop_CL7_per_fish_f20={};
fish=unique(idx_Fish_f20);
for clust=1:size(clustersF,1)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_f20_CL7_cleaned_cell.(clustersF{clust,1}));
for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_f20.(regionName).idx,temp_idx))/sum(ismember(PerBrainRegions_f20.(regionName).idx,intersect(idx_rsq_test_f20short_cleaned,temp_idx_fish)));
end

Prop_CL7_per_fish_f20{clust}(tempfish,:)=ratioperRegion;
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 
counter=1;
figure;
for clust=1:size(clustersF,1)

for i=1:length(fish)
   
   subplot(3,3,counter); plot(Prop_CL7_per_fish_f20{clust}(i,:),'.');set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
     
    hold on;
    %subplot(3,3,counter);bar(mean(Prop_CL7_per_fish_f20{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])

end
 
counter=counter+1;
end


counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);bar(mean(Prop_CL7_per_fish_f20{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end


%%% the box plot is very nice.
counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);boxplot(Prop_CL7_per_fish_f20{clust});set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end
%%
All_Prop_CL7_per_fish={};

for j=1:size(clustersF,1)
    All_Prop_CL7_per_fish.f20.(clustersF{j,1})=Prop_CL7_per_fish_f20{j};  
end

clearvars -except All_Prop_CL7_per_fish clustersF clustersS datasets RegionList

%%% the box plot is very nice.
counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);boxplot(All_Prop_CL7_per_fish.f20.(clustersF{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end


%%
%%% for s20

load('final_S20_step1.mat','idx_Fish_s20');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
BrainReg_S20=load('BrainReg_S20.mat');

PerBrainRegions_s20=BrainReg_S20.PerBrainRegions;

idx_rsq_test_s20short_cleaned=s20_cleaned_idxs.idx_rsq_test_s20short_cleaned;
idx_s20_multisense_cleaned=s20_cleaned_idxs.idx_multisense_cleaned;
idx_rsq_Mov_cleaned=s20_cleaned_idxs.idx_rsq_Mov_cleaned;
clust_s20_CL7_cleaned=s20_cleaned_idxs.clust_s20_CL7_cleaned;

clust_s20_CL7_cleaned_cell={};
clust=fieldnames(clust_s20_CL7_cleaned);
for j=1:size(clustersS,1)
 clust_s20_CL7_cleaned_cell.(clustersS{j,1})=clust_s20_CL7_cleaned.(clust{j});   
end    



%%% to look at the proportion per brain region to see if it works
%%% NOTE: i am basing the proportion on the response rsq filtered ROIs, not in
%%% the total number of ROIs in that brian region. it could be intersting
%%% to look at both.

counter=1;
figure;
for j=1:size(clustersS,1)
    
    idx_temp=clust_s20_CL7_cleaned_cell.(clustersS{j,1});
    
    for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_s20.(regionName).idx,idx_temp))/sum(ismember(PerBrainRegions_s20.(regionName).idx,idx_rsq_test_s20short_cleaned));
end
subplot(3,3,counter);bar(ratioperRegion);set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersS{j,1})); ylim([0 1])
    
   %hold on; 
   counter=counter+1;
end

%%% now i need to do it an average per fish

Prop_CL7_per_fish_s20={};
fish=unique(idx_Fish_s20);
for clust=1:size(clustersS,1)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s20==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_s20_CL7_cleaned_cell.(clustersS{clust,1}));
for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_s20.(regionName).idx,temp_idx))/sum(ismember(PerBrainRegions_s20.(regionName).idx,intersect(idx_rsq_test_s20short_cleaned,temp_idx_fish)));
end

Prop_CL7_per_fish_s20{clust}(tempfish,:)=ratioperRegion;
end
end

%%

for j=1:size(clustersS,1)
    All_Prop_CL7_per_fish.s20.(clustersS{j,1})=Prop_CL7_per_fish_s20{j};  
end

clearvars -except All_Prop_CL7_per_fish clustersS clustersS datasets RegionList

%%% the box plot is very nice.
counter=1; 
figure;
for clust=1:size(clustersS,1)
subplot(3,3,counter);boxplot(All_Prop_CL7_per_fish.s20.(clustersS{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersS{clust,1})); ylim([0 1])
 counter=counter+1;
end



%%
%%% to plot them in the same order

counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);boxplot(All_Prop_CL7_per_fish.f20.(clustersF{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end

counter=1; 
figure;
for clust=1:size(clustersF,1)
    
subplot(3,3,counter);boxplot(All_Prop_CL7_per_fish.s20.(clustersF{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end


%%
%%%% and together??

counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);
boxplot(All_Prop_CL7_per_fish.f20.(clustersF{clust}),'Colors','b');set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
hold on;
boxplot(All_Prop_CL7_per_fish.s20.(clustersF{clust}),'Colors','r');set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
counter=counter+1;
end



