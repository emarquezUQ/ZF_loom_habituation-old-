
%%%% this script is to look at the max responses per loom but in a
%%%% different way. I will check their distribution of strength and see in which
%%%% clusters the datasets differ the most


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


%% for f20
%%%for f20

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
%BrainReg_F20=load('BrainReg_F20.mat');

%PerBrainRegions_f20=BrainReg_F20.PerBrainRegions;

%%

%%% this is to get the ZS of the ROIs of each cluster per fish and per brain region of interest. 

clust_f20_CL4_cleaned=f20_cleaned_idxs.clust_f20_CL4_cleaned;

clust_f20_CL4_cleaned_cell={};
clust=fieldnames(clust_f20_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL4_cleaned_cell.(clustersF{j,1})=clust_f20_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_f20=struct;
ZS_CL4_per_fishNbrain_f20=struct;
fish=unique(idx_Fish_f20);
edges=[0:0.25:15];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f20.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
        for k=1:30

        temp_max_resp{1,k}= (max(ZS_f20(temp_idx,loom_moments{1,k})'))-(ZS_f20(temp_idx,Loomf20_onset_idx(k)+1))'; 
        
        [count_temp_max_resp{1,k},~]=histcounts(temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_f20.maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=temp_max_resp;
ZS_CL4_per_fishNbrain_f20.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_f20.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_f20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

end
end
end

%%% to check it worked
 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;


for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
    
for i=[1 2 3 11]
    
    subplot(4,4,counter);
   plot(ZS_CL4_per_fishNbrain_f20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
   counter=counter+1; 
end


end
sgtitle(RegionList{brain});

end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS ZS_CL4_per_fishNbrain_f20


%% for f60
%%%for f60

load('final_F60_step1_2.mat','ZS_f60','idx_Fish_f60','ZS_short_F60');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');


%%

%%% this is to get the ZS of the ROIs of each cluster per fish and per brain region of interest. 

clust_f60_CL4_cleaned=f60_cleaned_idxs.clust_f60_CL4_cleaned;

clust_f60_CL4_cleaned_cell={};
clust=fieldnames(clust_f60_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f60_CL4_cleaned_cell.(clustersF{j,1})=clust_f60_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_f60=struct;
ZS_CL4_per_fishNbrain_f60=struct;
fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...

edges=[0:0.25:15];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f60_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f60.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
        for k=1:30

        temp_max_resp{1,k}= (max(ZS_f60(temp_idx,ZS_short_F60(loom_moments{1,k}))'))-(ZS_f60(temp_idx,ZS_short_F60(Loomf20_onset_idx(k)+1)))'; 
        
        [count_temp_max_resp{1,k},~]=histcounts(temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_f60.maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=temp_max_resp;
ZS_CL4_per_fishNbrain_f60.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_f60.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_f60.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

end
end
end

%%% to check it worked
 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;


for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
    
for i=[1 2 3 11]
    
    subplot(4,4,counter);
   plot(ZS_CL4_per_fishNbrain_f60.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
   counter=counter+1; 
end


end
sgtitle(RegionList{brain});

end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS ZS_CL4_per_fishNbrain_f20 ZS_CL4_per_fishNbrain_f60


%% for s20

%%% now for s20

load('final_S20_step1.mat','ZS_s20','idx_Fish_s20');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');

%%

%%% this is to get the ZS of the ROIs of each cluster per fish and per brain region of interest. 

clust_s20_CL4_cleaned=s20_cleaned_idxs.clust_s20_CL4_cleaned;

clust_s20_CL4_cleaned_cell={};
clust=fieldnames(clust_s20_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s20_CL4_cleaned_cell.(clustersS{j,1})=clust_s20_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_s20=struct;
ZS_CL4_per_fishNbrain_s20=struct;
fish=unique(idx_Fish_s20);

edges=[0:0.25:15];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s20_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s20.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
        for k=1:30

        temp_max_resp{1,k}= (max(ZS_s20(temp_idx,S_trim(loom_moments{1,k}))'))-(ZS_s20(temp_idx,S_trim(Loomf20_onset_idx(k)+1)))'; 
        
        [count_temp_max_resp{1,k},~]=histcounts(temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_s20.maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=temp_max_resp;
ZS_CL4_per_fishNbrain_s20.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_s20.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_s20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

end
end
end

%%% to check it worked
 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;


for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
    
for i=[1 2 3 11]
    
    subplot(4,4,counter);
   plot(ZS_CL4_per_fishNbrain_s20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
   counter=counter+1; 
end


end
sgtitle(RegionList{brain});

end

clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS ZS_CL4_per_fishNbrain_f20 ZS_CL4_per_fishNbrain_f60 ZS_CL4_per_fishNbrain_s20


%% for s60


%%% now for s60

load('final_S60_step1.mat','ZS_s60','idx_Fish_s60','ZS_short_S60');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');




%%

%%% this is to get the ZS of the ROIs of each cluster per fish and per brain region of interest. 

clust_s60_CL4_cleaned=s60_cleaned_idxs.clust_s60_CL4_cleaned;

clust_s60_CL4_cleaned_cell={};
clust=fieldnames(clust_s60_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s60_CL4_cleaned_cell.(clustersS{j,1})=clust_s60_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_s60=struct;
ZS_CL4_per_fishNbrain_s60=struct;
fish=unique(idx_Fish_s60);

edges=[0:0.25:15];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s60_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s60.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN(1,30)}; 
    
    %count_temp_max_resp={NaN(1,60)};
    
else
    
        for k=1:30

        temp_max_resp{1,k}= (max(ZS_s60(temp_idx,ZS_short_S60(S_trim(loom_moments{1,k})))'))-(ZS_s60(temp_idx,ZS_short_S60(S_trim(Loomf20_onset_idx(k)+1))))'; 
        
        [count_temp_max_resp{1,k},~]=histcounts(temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_s60.maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=temp_max_resp;
ZS_CL4_per_fishNbrain_s60.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_s60.countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_s60.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

end
end
end

%%% to check it worked
 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;


for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for i=[1 2 3 11]
    
    subplot(4,4,counter);
   plot(ZS_CL4_per_fishNbrain_s60.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
   counter=counter+1; 
end


end
sgtitle(RegionList{brain});

end

%% plotting together


for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;


for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for i=[1 2 3 11]
    
    subplot(4,4,counter);
   plot(ZS_CL4_per_fishNbrain_f20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold on;
   plot(ZS_CL4_per_fishNbrain_f60.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold on;
   plot(ZS_CL4_per_fishNbrain_s20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold on;
   plot(ZS_CL4_per_fishNbrain_s60.MeancountMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold off;
   
   
   counter=counter+1; 
end


end
sgtitle(RegionList{brain});

end

clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments S_trim RegionList PerBrainRegions datasets clustersF clustersS ZS_CL4_per_fishNbrain_f20 ZS_CL4_per_fishNbrain_f60 ZS_CL4_per_fishNbrain_s20 ZS_CL4_per_fishNbrain_s60


save('Max_response_perBrain_all_CL4_corrected_distribution.mat');

%%% it seems that the only interesting part is in the tectum... the other
%%% ones look a bit too noisy to really say something. something weird is
%%% that I aslo see differences in the first loom when it should be the
%%% same for looms of the same speed... although as this is rsq filtered
%%% ROIs there could be a difference there... Maybe I need to checkit with
%%% just filtering with the df/f instead. 
