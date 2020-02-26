
%%%% this script is to get the differences in loom responses across
%%%% datasets and cluster types. I will take the peak of each loom and compare it to the
%%%% first loom response in each ROI for all 4 datasets. Then I will get
%%%% the mean of each fish for each dataset. 





%%% first with f20


%%%% first you need to load what you need

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

load('BrainReg_F20.mat')
load('final_F20_step1.mat','MatFiles_f20','ZS_f20','idx_Fish_f20','idx_Plane_f20','rawregressF20','idx_rsq_test_f20short','High_corr_Nb_f20','High_corr_Nb_f20_short');
load('Zbrain_Masks.mat');
load('final_F20_step1.mat','gooodmaps');
load('f20_cleaned_idxs.mat')
load('allfish_f20looms_tailmov.mat');

%%
 %%% this is to get the onset and the loom moments

%%% to take the loom onsets based on the movie

%%% for f20
Loomf20_onset=zeros(6720,1);
startpoints=[30,254,478]; %%in seconds
loom_times=[0,22,40,60,78,100,120,140,158,180]; %%in seconds
for p=1:3
for k=1:10
    Loomf20_onset(startpoints(p)*10+loom_times(k)*10)=2;
end
end

figure; plot(Loomf20_onset);
%%

Loomf20_onset_idx=find(Loomf20_onset==2);
Loomf20_onset_idx=round(Loomf20_onset_idx/5);
Loomf20_onset_idx(Loomf20_onset_idx<1)=[];
Loomf20_onset_idx(Loomf20_onset_idx>(length(allfish_f20looms_tailmov)-1))=[];        
  
%%
%%% this is to try to get the idx to normalize with each loom. 


loom_moments={};
for i=1:30
    if i==1||i==11||i==21
    loom_moments{i}=(Loomf20_onset_idx(i)-59):Loomf20_onset_idx(i+1);
    elseif i==10||i==20||i==30
      loom_moments{i}=Loomf20_onset_idx(i)+1:Loomf20_onset_idx(i)+28;  
    else
     loom_moments{i}=Loomf20_onset_idx(i)+1:Loomf20_onset_idx(i+1);   
    end
end

%%% to check
for i=1:30
temp(1,i)=length(loom_moments{i});

end
sum(temp(1,:))

%%

%%% here i was trying to normalize each ROI based on the response to the
%%% first loom but it doesnt work well for all the clusters. for example it
%%% doesnt work with the nonhab cause some of the spikes afterwards are
%%% stronger that the first one (maybe there was a small drift in that
%%% ROI). the results are pretty bad in those cases. it does work for
%%% fasthab cluster though. 


Max_resp_f20=[];

for i=1:length(idx_rsq_test_f20short_cleaned)

for k=1:30

   temp_Max_resp(1,k)= max(ZS_f20(idx_rsq_test_f20short_cleaned(i),loom_moments{1,k}))-min(ZS_f20(idx_rsq_test_f20short_cleaned(i),Loomf20_onset_idx(k)+1));
end


%Max_resp_f20(i,:)=temp_Max_resp(i,:)/temp_Max_resp(i,1);
Max_resp_f20(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);

end

figure;plot(mean(Max_resp_f20));

figure;imagesc(Max_resp_f20);


temp_idx=find(ismember(idx_rsq_test_f20short_cleaned,clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned));




figure;plot(mean(Max_resp_f20(temp_idx,:)));

figure;imagesc(Max_resp_f20(temp_idx,:));


figure;imagesc(ZS_f20(clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned,:));
figure;imagesc(ZS_f20(idx_rsq_test_f20short_cleaned(temp_idx),:));


figure;
for i=1:length(Max_resp_f20(temp_idx))
   plot(Max_resp_f20(temp_idx(i),:));
hold on;
    
end


%%

%%% NOTE... it seems that i cannot do it to every individual ROI. cause in
%%% some of the clusters the variability and strength makes it harder to
%%% normalize. so I will try per fish.

%%% this is how the responses of a single fish varies in one cluster to
%%% show what I am mentioning above.

temp_idx2=intersect(clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned,find(idx_Fish_f20==1));

figure;plot(mean(ZS_f20(temp_idx2,:)));

figure;
for i=1:length(temp_idx2)
plot(ZS_f20(temp_idx2(i),:));
hold on;
end

figure;imagesc(ZS_f20(temp_idx2,:));




%%

%%% trying first to normalize each ROI with in each cluster. I am
%%% normalizing to the highest value not to the first loom. 

%%% the result is interesting...

norm_rsq=[];
for i=1:length(idx_rsq_test_f20short_cleaned)
    for k=1:30
temp_norm_rsq(1,k)=max(ZS_f20(idx_rsq_test_f20short_cleaned(i),loom_moments{1,k}))-min(ZS_f20(idx_rsq_test_f20short_cleaned(i),Loomf20_onset_idx(k)+1));

    end
norm_rsq(i,:)=temp_norm_rsq(1,:)/max(temp_norm_rsq(1,:));
end


figure;plot(mean(norm_rsq));

figure;imagesc(norm_rsq);




temp_idx=find(ismember(idx_rsq_test_f20short_cleaned,clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned));


figure;plot(mean(norm_rsq(temp_idx,:)));

figure;imagesc(norm_rsq(temp_idx,:));

figure;plot(mean(ZS_f20(clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned,:)));
figure;imagesc(ZS_f20(clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned,:), [-0.5 4]);colormap hot
figure;imagesc(ZS_f20(idx_rsq_test_f20short_cleaned(temp_idx),:), [-0.5 4]);colormap hot





%%
%%% testing some things

%%% what if i normalize if ROI before?
%%% although it might not be necessary if I then take the mean of each
%%% fish. 


%%% note: the normalize function only works in matlab2018...

figure;imagesc(ZS_f20(idx_rsq_test_f20short_cleaned,:));



figure;imagesc(normalize(ZS_f20(idx_rsq_test_f20short_cleaned,:)));


figure;plot(normalize(ZS_f20(idx_rsq_test_f20short_cleaned(i),:),'range'));


normZS_rsq=[];
temp_norm_rsq=[];
for i=1:length(idx_rsq_test_f20short_cleaned)
temp_norm_rsq(1,:)=normalize(ZS_f20(idx_rsq_test_f20short_cleaned(i),:),'range');
normZS_rsq(i,:)=temp_norm_rsq;
end


figure;imagesc(normZS_rsq);


%%% checking with the nonhab cluster
figure;imagesc(normZS_rsq(temp_idx,:));



figure;plot(mean(normZS_rsq(temp_idx,:)));

figure;
for i=1:length(temp_idx)
plot(normZS_rsq(temp_idx(i),:));
hold on;
end

%%
%%% trying again to normalize each ROI with in each cluster

%%% the result is ... not working yet

norm_rsq2=[];
for i=1:length(normZS_rsq)
    for k=1:30
temp_norm_rsq(1,k)=max(normZS_rsq(i,loom_moments{1,k}))-min(normZS_rsq(i,Loomf20_onset_idx(k)+1));

    end
norm_rsq2(i,:)=temp_norm_rsq(1,:)/max(temp_norm_rsq(1,:));
end


figure;plot(mean(norm_rsq2));

figure;imagesc(norm_rsq2);




temp_idx=find(ismember(idx_rsq_test_f20short_cleaned,clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned));


figure;plot(mean(norm_rsq2(temp_idx,:)));

figure;imagesc(norm_rsq2(temp_idx,:));

%%




%%% this is to get the means of each cluster per fish. 


clust_f20_CL7_cleaned2=struct2cell(clust_f20_CL7_cleaned);

Means_CL7_per_fish_f20={}
fish=unique(idx_Fish_f20);
for clust=1:length(clust_f20_CL7_cleaned2)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_f20_CL7_cleaned2{clust,1});
temp_mean=mean(ZS_f20(temp_idx,:));

Means_CL7_per_fish_f20{clust}(tempfish,:)=temp_mean;
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_f20_CL7_cleaned2)
figure;
for i=1:length(fish)
   plot(Means_CL7_per_fish_f20{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Means_CL7_per_fish_f20{1,clust}))

end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f20_perfish={};

for clust=1:length(clust_f20_CL7_cleaned2)

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL7_per_fish_f20{1,clust}(i,loom_moments{1,k}))-min(Means_CL7_per_fish_f20{1,clust}(i,Loomf20_onset_idx(k)+1));
end


Max_resp_f20_perfish{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_f20_CL7_cleaned2)
figure;
for i=1:length(fish)
   plot(Max_resp_f20_perfish{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Max_resp_f20_perfish{1,clust}))

end



%%

%%% to trim the slow loom movies

S_trim=[1:448 453:901 906:1352];


%%


%%% this is to compare the datasets


datasets=['f20'; 'f60'; 's20'; 's60'];
clusters=fieldnames(mean_CL7_f20);





