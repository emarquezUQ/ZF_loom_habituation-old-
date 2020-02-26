

%%%% this script is to look at the max responses per loom but in a
%%%% different way. I will check their distribution of strength and see in which
%%%% clusters the datasets differ the most


%%% it seems that I was getting an artifact... the resposnes to the first
%%% loom also show differences between different ISI but same speed. this
%%% response should be very similar because the first loom is the same and
%%% apears at the same time with in a speed group. I think it is becasue
%%% the Z scoring is different for the longer movies. so I need to solve
%%% this... I think I need to normalize them based on the max response of
%%% the first loom. so I made this script to try to correct for that. it is
%%% based in Max_resposne_perBrain_CL4_corrected_distribution.m


%%%% but the results look weird... I am not convinced it worked. 


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
edges=[0:0.01666:1];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f20.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
    k=1;
        temp_max_resp1= (max(ZS_f20(temp_idx,loom_moments{1,k})'))-(ZS_f20(temp_idx,Loomf20_onset_idx(k)+1))';
        temp_min=min(temp_max_resp1);
        temp_max=max(temp_max_resp1);
        
        temp_max_resp{1,k}=(max(ZS_f20(temp_idx,loom_moments{1,k})'))-(ZS_f20(temp_idx,Loomf20_onset_idx(k)+1))';
        
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
    
        for k=2:30

        temp_max_resp{1,k}= (max(ZS_f20(temp_idx,loom_moments{1,k})'))-(ZS_f20(temp_idx,Loomf20_onset_idx(k)+1))'; 
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_f20.Norm_maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_temp_max_resp;
ZS_CL4_per_fishNbrain_f20.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_f20.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_f20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

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
   plot(ZS_CL4_per_fishNbrain_f20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
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

edges=[0:0.01666:1];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f60_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f60.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
    k=1;
        temp_max_resp1= (max(ZS_f60(temp_idx,loom_moments{1,k})'))-(ZS_f60(temp_idx,ZS_short_F60(Loomf20_onset_idx(k)+1)))';
        temp_min=min(temp_max_resp1);
        temp_max=max(temp_max_resp1);
        
        temp_max_resp{1,k}=(max(ZS_f60(temp_idx,loom_moments{1,k})'))-(ZS_f60(temp_idx,ZS_short_F60(Loomf20_onset_idx(k)+1)))';
        
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
    
        for k=2:30

        temp_max_resp{1,k}= (max(ZS_f60(temp_idx,loom_moments{1,k})'))-(ZS_f60(temp_idx,ZS_short_F60(Loomf20_onset_idx(k)+1)))'; 
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_f60.Norm_maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_temp_max_resp;
ZS_CL4_per_fishNbrain_f60.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_f60.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_f60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

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
   plot(ZS_CL4_per_fishNbrain_f60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
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

edges=[0:0.01666:1];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s20_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s20.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
    k=1;
        temp_max_resp1= (max(ZS_s20(temp_idx,loom_moments{1,k})'))-(ZS_s20(temp_idx,S_trim(Loomf20_onset_idx(k)+1)))';
        temp_min=min(temp_max_resp1);
        temp_max=max(temp_max_resp1);
        
        temp_max_resp{1,k}=(max(ZS_s20(temp_idx,loom_moments{1,k})'))-(ZS_s20(temp_idx,S_trim(Loomf20_onset_idx(k)+1)))';
        
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
    
        for k=2:30

        temp_max_resp{1,k}= (max(ZS_s20(temp_idx,loom_moments{1,k})'))-(ZS_s20(temp_idx,S_trim(Loomf20_onset_idx(k)+1)))'; 
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_s20.Norm_maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_temp_max_resp;
ZS_CL4_per_fishNbrain_s20.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_s20.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_s20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

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
   plot(ZS_CL4_per_fishNbrain_s20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
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

edges=[0:0.01666:1];

for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s60_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s60.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
    k=1;
        temp_max_resp1= (max(ZS_s60(temp_idx,loom_moments{1,k})'))-(ZS_s60(temp_idx,ZS_short_S60(S_trim(Loomf20_onset_idx(k)+1))))';
        temp_min=min(temp_max_resp1);
        temp_max=max(temp_max_resp1);
        
        temp_max_resp{1,k}=(max(ZS_s60(temp_idx,loom_moments{1,k})'))-(ZS_s60(temp_idx,ZS_short_S60(S_trim(Loomf20_onset_idx(k)+1))))';
        
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
    
        for k=2:30

        temp_max_resp{1,k}= (max(ZS_s60(temp_idx,loom_moments{1,k})'))-(ZS_s60(temp_idx,ZS_short_S60(S_trim(Loomf20_onset_idx(k)+1))))'; 
        norm_temp_max_resp{1,k}=(temp_max_resp{1,k}-temp_min)/(temp_max-temp_min);
        
        [norm_count_temp_max_resp{1,k},~]=histcounts(norm_temp_max_resp{1,k},edges);
        end


end

ZS_CL4_per_fishNbrain_s60.Norm_maxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_temp_max_resp;
ZS_CL4_per_fishNbrain_s60.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}=norm_count_temp_max_resp;



end
end
end



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
 for k=1:30
     
     temp_count=zeros(length(fish),60);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_CL4_per_fishNbrain_s60.Norm_countMaxresponsePerloom.(RegionList{brain}).(clustersF{clust}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_CL4_per_fishNbrain_s60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,k}=temp_mean; 

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
   plot(ZS_CL4_per_fishNbrain_s60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   
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
    
   plot(ZS_CL4_per_fishNbrain_f20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold on;
   plot(ZS_CL4_per_fishNbrain_f60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold on;
   plot(ZS_CL4_per_fishNbrain_s20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold on;
   plot(ZS_CL4_per_fishNbrain_s60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold off;
   
   
   counter=counter+1; 
end


end
sgtitle(RegionList{brain});

end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments S_trim RegionList PerBrainRegions datasets clustersF clustersS ZS_CL4_per_fishNbrain_f20 ZS_CL4_per_fishNbrain_f60 ZS_CL4_per_fishNbrain_s20 ZS_CL4_per_fishNbrain_s60


%save('Max_response_perBrain_all_CL4_corrected_distribution2.mat');

%%% it seems that the only interesting part is in the tectum... the other
%%% ones look a bit too noisy to really say something. something weird is
%%% that I aslo see differences in the first loom when it should be the
%%% same for looms of the same speed... although as this is rsq filtered
%%% ROIs there could be a difference there... Maybe I need to checkit with
%%% just filtering with the df/f instead. 


%% normalizing based on the area under the curve



edges=[0:0.01667:1];

for brain=1:length(RegionList)



for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for i=1:30
    
   
   temp_mean=ZS_CL4_per_fishNbrain_f20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}; 
   temp_area=trapz(edges,temp_mean);
   temp_mean_norm=temp_mean./temp_area;
   ZS_CL4_per_fishNbrain_f20.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}=temp_mean_norm;
   
   temp_mean=ZS_CL4_per_fishNbrain_f60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}; 
   temp_area=trapz(edges,temp_mean);
   temp_mean_norm=temp_mean./temp_area;
   ZS_CL4_per_fishNbrain_f60.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}=temp_mean_norm;
   
   temp_mean=ZS_CL4_per_fishNbrain_s20.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}; 
   temp_area=trapz(edges,temp_mean);
   temp_mean_norm=temp_mean./temp_area;
   ZS_CL4_per_fishNbrain_s20.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}=temp_mean_norm;
   
   temp_mean=ZS_CL4_per_fishNbrain_s60.Mean_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}; 
   temp_area=trapz(edges,temp_mean);
   temp_mean_norm=temp_mean./temp_area;
   ZS_CL4_per_fishNbrain_s60.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}=temp_mean_norm; 
   
   trapz(edges, temp_mean_norm)
    
end


end


end

%% to plot them together after being normalized. 


edges=[0:0.01667:1];
counter=1; 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for i=[1 2 3 4 5 11]
    
   subplot(4,6,counter);
    
   plot(ZS_CL4_per_fishNbrain_f20.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i});
   hold on;
   
   plot(ZS_CL4_per_fishNbrain_f60.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i});
   hold on;
   
   plot(ZS_CL4_per_fishNbrain_s20.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i});
   hold on;
   
   plot(ZS_CL4_per_fishNbrain_s60.MeanNorm_norm_countMaxresponsePerloomPerfish.(RegionList{brain}).(clustersF{clust}){1,i}); 
   hold on;
   
   counter=counter+1;  
   
end
end
sgtitle(RegionList{brain});
end


