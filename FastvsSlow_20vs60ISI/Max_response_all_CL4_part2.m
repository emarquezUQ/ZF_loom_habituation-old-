%%% this script is a follow up of the Max_response_all_CL4.m and is just to
%%% add a timepoint 0 for my max response graphs. I am taking the onset + 1
%%% of the first loom. 

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab


load('All_means_maxResp_CL4_normalized.mat')


%%% for f20

load('final_F20_step1.mat','idx_Fish_f20');

fish=unique(idx_Fish_f20);


Max_resp_f20_perfish2={};

for clust=1:4

for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fish_f20{1,clust}(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fish_f20{1,clust}(i,Loomf20_onset_idx(1)+1));
    
for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fish_f20{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_f20{1,clust}(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint zero
end


Max_resp_f20_perfish2{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2); %%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
end
end



%%

%%% for f60

load('final_F60_step1_2.mat','idx_Fish_f60');

fish=unique(idx_Fish_f60);

fish(find(fish==47))=[]; %%% cause I also took out fish 47...


Max_resp_f60_perfish2={};

for clust=1:4

for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fish_f60{1,clust}(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fish_f60{1,clust}(i,Loomf20_onset_idx(1)+1));
    
for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fish_f60{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_f60{1,clust}(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint
end


Max_resp_f60_perfish2{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2); %%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
end
end

%%
%%% for s20

load('final_S20_step1.mat','idx_Fish_s20');

fish=unique(idx_Fish_s20);

Max_resp_s20_perfish2={};

for clust=1:4

for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fish_s20{1,clust}(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fish_s20{1,clust}(i,Loomf20_onset_idx(1)+1));
    
for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fish_s20{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_s20{1,clust}(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint
end


Max_resp_s20_perfish2{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2); %%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
end
end

%%

%%% now for s60

load('final_S60_step1.mat','idx_Fish_s60');
fish=unique(idx_Fish_s60);

Max_resp_s60_perfish2={};

for clust=1:4

for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fish_s60{1,clust}(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fish_s60{1,clust}(i,Loomf20_onset_idx(1)+1));
    
for k=1:30

   temp_Max_resp(1,k+1)= max(Means_CL4_per_fish_s60{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_s60{1,clust}(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint
end


Max_resp_s60_perfish2{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2); %%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
end
end


%%
%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:4
figure;
for i=1:length(fish)
   plot(Max_resp_s60_perfish2{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Max_resp_s60_perfish2{1,clust}))

end


clearvars -except Loomf20_onset loom_moments Loomf20_onset_idx S_trim Means_CL4_per_fish_f20 Max_resp_f20_perfish2 Means_CL4_per_fish_f60 Max_resp_f60_perfish2 Means_CL4_per_fish_s20 Max_resp_s20_perfish2 Means_CL4_per_fish_s60 Max_resp_s60_perfish2

%%
%%% this is to compare the datasets


datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);


%%

%%% organizing the data
All_means_maxResp_CL4_normalized2={};


counter=1;
figure;
for data=1:size(datasets,1)
  datasetname=strcat('Max_resp_',datasets(data,:),'_perfish2');
  tempVar=eval(datasetname);

if ~isempty(regexp(datasets(data,:),'f'))
      
    for j=1:size(clustersF,1)


    All_means_maxResp_CL4_normalized2.(datasets(data,:)).(clustersF{j,1})=nanmean(tempVar{1,j});

    subplot(4,7,counter);plot(All_means_maxResp_CL4_normalized2.(datasets(data,:)).(clustersF{j,1}));

    counter=counter+1
    end

else

    for j=1:size(clustersS,1)


    All_means_maxResp_CL4_normalized2.(datasets(data,:)).(clustersS{j,1})=nanmean(tempVar{1,j});

    subplot(4,7,counter);plot(All_means_maxResp_CL4_normalized2.(datasets(data,:)).(clustersS{j,1}));

    counter=counter+1

 end

    
end

end


clear tempVar  datasetname

All_means_maxResp_CL4_normalized2=struct2cell(All_means_maxResp_CL4_normalized2);


for j=1:4
figure;
for i=1:4
   plot(All_means_maxResp_CL4_normalized2{i,1}.(clustersF{j,1}))

hold on;
end 



end
        
        
 clearvars -except Loomf20_onset loom_moments Loomf20_onset_idx S_trim Means_CL4_per_fish_f20 Max_resp_f20_perfish2 Means_CL4_per_fish_f60 Max_resp_f60_perfish2 Means_CL4_per_fish_s20 Max_resp_s20_perfish2 Means_CL4_per_fish_s60 Max_resp_s60_perfish2 All_means_maxResp_CL4_normalized2
       
save('All_means_maxResp_CL4_normalized2.mat');

