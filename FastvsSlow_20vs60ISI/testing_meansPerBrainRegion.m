
%%%% this script is to get the means per brain regions

%%% for f20: I will try  with the fasthab and nonhab clusters.  

%%% only for some of the brain areas: Pallium, Thalamus, Tectum, hindbrain



%%%%% the results are not that interesting... not much different that I
%%%%% already have with figure 3 cause there is not that much variability
%%%%% in the nonhab cluster so there are more differences in fasthab.
%%%%% Maybe a proportion figure for fig 3 would be better

load('f20_cleaned_idxs.mat');

load('All_More_BrainReg2.mat');

load('final_F20_step1.mat','idx_Fish_f20');

load('means_F20_CL4n7.mat','mean_CL4_f20');

load('Max_response_perBrain_all_CL4_corrected.mat','Means_CL4_per_fishNbrain_f20');


%RegionList={'Pallium','Thalamus','Pretectum','Tectum','Hindbrain'};

%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};


clustersF=fieldnames(mean_CL4_f20);

clust_f20_CL4_cleaned_cell={};
clust=fieldnames(clust_f20_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL4_cleaned_cell.(clustersF{j,1})=clust_f20_CL4_cleaned.(clust{j});   
end  


fish=unique(idx_Fish_f20);


%%% for fasthab 

RegionList={'Pallium','Thalamus','Tectum','Tegmentum','Hindbrain'};

figure('Position',[100 0 900 900]);
counter=1;
for brain=1:length(RegionList)


for clust=2%length(fieldnames(clust_f20_CL4_cleaned_cell))  %%% i changed it to 2 to see only fasthab
subplot(5,1,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f20{1,clust}))

end
end


%%


%%% for nonhab


RegionList={'Pallium','Habenula','Thalamus','Pretectum','Tectum'};
figure('Position',[100 0 900 900]);
counter=1;
for brain=1:length(RegionList)


for clust=1%length(fieldnames(clust_f20_CL4_cleaned_cell))  %%% i changed it to 1 to see only nonhab
subplot(5,1,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f20{1,clust}))

end
end