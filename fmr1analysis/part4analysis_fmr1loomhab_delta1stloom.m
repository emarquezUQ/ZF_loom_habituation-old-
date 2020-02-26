
%%% this script is to look at the amount of ROIs that respond to the first loom.    
%%% the idea is to check if the reason why i get less ROIs in the fmr1 and het fish     
%%% is due to a filtering with the r2 rather than a lower number of respondent neurons  

%%% I will select this ROIs with dF/F to a 2SD response and maybe a simpble
%%% correlation. 

%%% I will select the df/f threshold from the wildtype fish. 



%%%% then there is a second part where i reclustered the ROIs of each
%%%% genotype and did a hierarchical clustering based on a correlation to
%%%% merge the ones that are similar. It was a bit interesting but I dont
%%%% think I go new information. 

load('s20_fmr1_loomhab_CN.mat','MatFiles','ZS_CN');


load('s20_fmr1_loomhab_CN_post_rsq01.mat','idx_rsq1');

%%% to get the clusters from the Kmeans done to ALL the ROIs
load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN');

load('s20_good_idx_Fish.mat','idx_Fish');
idx_Fish_cat=categorical(idx_Fish);

load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');


load('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','idx_rsq');

load('s20_fmr1_loomhab_CN_part2_High_corr_Nb.mat','High_corr_Nb');

 %%%% i also need to generate the list of fish saved in the
 %%%% fmr1loomhab_lists.m file
%%% or load idx_temp2=fmr1; idx_temp4=controls and idx_temp5=hets
load('s20_fmr1_loomhab_CN_part3.mat','idx_temp2','idx_temp4','idx_temp5'); 


%%
%%% first getting the ROIs of each genotype

idx_het1=ismember(idx_Fish,list1);
idx_het1=find(idx_het1);
idx_fmr1=ismember(idx_Fish,list2);
idx_fmr1=find(idx_fmr1);
idx_het2=ismember(idx_Fish,list3);
idx_het2=find(idx_het2);
idx_wt=ismember(idx_Fish,list4);
idx_wt=find(idx_wt);
idx_hets=union(idx_het1,idx_het2);


%% getting the delta of the first loom

%%% this is if I want to get first the ROIs that reacted when the first
%%% loom presentation happen based on the strenght of the response
%%% filtering by 3 units of the z-score. then I get the delta of the first
%%% response. all this in wild type fish
idx_resp_wt=[];
delta_1_wt=[];
for i=1:size(ZS_CN(idx_wt),1)

%if  max(ZS_CN(i,60:80))>2   
idx_resp_wt(i)=max(ZS_CN(idx_wt(i),60:80))>3;
%elseif min(ZS_CN(i,60:80))<-2 
%idx_resp(i)=min(ZS_CN(i,60:80))<-2;
%end

delta_1_wt(i)=max(ZS_CN(idx_wt(i),60:80))-mean(ZS_CN(idx_wt(i),10:59));
% delta_2(i)=max(ZS_CN(i,420:440))-mean(ZS_CN(i,10:59));
% deltahab(i)=delta_1(i)-delta_2(i);

  
end

idx_resp_wt=idx_wt(find(idx_resp_wt));

figure;imagesc(ZS_CN(idx_resp_wt,:)); colormap('hot');
figure;histogram(delta_1_wt);

 
%%% then I get the SD of the difference between 1st and 10th to select the 
%%% ROIs that had above or bellow a 2 SD difference

std_delta_1_wt=std(delta_1_wt);

idx_delta_pos=find(delta_1_wt>std_delta_1_wt*2);
idx_delta_neg=find(delta_1_wt<-std_delta_1_wt*2);
idx_delta=union(idx_delta_neg,idx_delta_pos);

histogram(delta_1_wt(idx_delta));

%%% there are no negative ones...

figure;imagesc(ZS_CN(idx_wt(idx_delta_pos),:));  colormap('hot');
figure;plot(mean(ZS_CN(idx_wt(idx_delta_pos),:)));  
  

%%% so i will only use the positive ones

deltahab_2std=delta_1_wt(idx_delta_pos);
Qs = quantile(deltahab_2std,[0.025 0.25 0.50 0.75 0.975]);


idx_delta_pos=idx_wt(idx_delta_pos); %%% to put it in the right indexing format

figure;imagesc(ZS_CN(idx_delta_pos,:));  colormap('hot');
figure;plot(mean(ZS_CN(idx_delta_pos,:)));

%%% this is for when I get the ANTs results to be able to plot the ROIs in
%%% Zbrain

% figure;
% scatter(ROI_temp2.f20(idx_delta_pos,1),ROI_temp2.f20(idx_delta_pos,2),10,deltahab_2std,'filled');colormap('jet');colorbar; caxis([Qs(1) Qs(5)]);
% figure;
% scatter3(ROI_temp2.f20(idx_delta_pos,1),ROI_temp2.f20(idx_delta_pos,2),ROI_temp2.f20(idx_delta_pos,3),10,deltahab_2std,'filled');colormap('jet');colorbar; caxis([Qs(1) Qs(5)]);
%  
%%%%% it looks ... 


%% now checking fmr1 and hets



idx_resp_fmr1=[];
delta_1_fmr1=[];
for i=1:size(ZS_CN(idx_fmr1),1)

%if  max(ZS_CN(i,60:80))>2   
idx_resp_fmr1(i)=max(ZS_CN(idx_fmr1(i),60:80))>3;
%elseif min(ZS_CN(i,60:80))<-2 
%idx_resp(i)=min(ZS_CN(i,60:80))<-2;
%end

delta_1_fmr1(i)=max(ZS_CN(idx_fmr1(i),60:80))-mean(ZS_CN(idx_fmr1(i),10:59));
% delta_2(i)=max(ZS_CN(i,420:440))-mean(ZS_CN(i,10:59));
% deltahab(i)=delta_1(i)-delta_2(i);

  
end

idx_resp_fmr1=idx_fmr1(find(idx_resp_fmr1));

figure;imagesc(ZS_CN(idx_resp_fmr1,:)); colormap('hot');
figure;histogram(delta_1_fmr1);

 
%%% then I  select the ROIs that had above or bellow a 2 SD difference of
%%% the delta

std_delta_1_fmr1=std(delta_1_fmr1); %%% the standard deviation is simiar to wt... 

idx_delta_pos_fmr1=find(delta_1_fmr1>std_delta_1_wt*2);  %%% i am using as a filter 2SD of the wild type delta
idx_delta_neg_fmr1=find(delta_1_fmr1<-std_delta_1_wt*2);
idx_delta_fmr1=union(idx_delta_neg_fmr1,idx_delta_pos_fmr1);

histogram(delta_1_fmr1(idx_delta_fmr1));

%%% there are no negative ones...

figure;imagesc(ZS_CN(idx_fmr1(idx_delta_pos_fmr1),:));  colormap('hot');
figure;plot(mean(ZS_CN(idx_fmr1(idx_delta_pos_fmr1),:)));  
  

%%% so i will only use the positive ones

deltahab_2std_fmr1=delta_1_fmr1(idx_delta_pos_fmr1);
Qs = quantile(deltahab_2std_fmr1,[0.025 0.25 0.50 0.75 0.975]);


idx_delta_pos_fmr1=idx_fmr1(idx_delta_pos_fmr1); %%% to put it in the right indexing format

figure;imagesc(ZS_CN(idx_delta_pos_fmr1,:));  colormap('hot');
figure;plot(mean(ZS_CN(idx_delta_pos_fmr1,:)));


%%% it seems that even like this wild types have a higher number of ROIs...
%%% what if check with the SD of the fmr1? it increases a bit the number of
%%% ROIs of fmr1 (from 21K to 24K aprox). but still bellow the 29K of wild
%%% types. this is not yet normalized per fish yet though. 

%%
%%% now lets check at the hets. 

idx_resp_hets=[];
delta_1_hets=[];
for i=1:size(ZS_CN(idx_hets),1)

%if  max(ZS_CN(i,60:80))>2   
idx_resp_hets(i)=max(ZS_CN(idx_hets(i),60:80))>3;
%elseif min(ZS_CN(i,60:80))<-2 
%idx_resp(i)=min(ZS_CN(i,60:80))<-2;
%end

delta_1_hets(i)=max(ZS_CN(idx_hets(i),60:80))-mean(ZS_CN(idx_hets(i),10:59));
% delta_2(i)=max(ZS_CN(i,420:440))-mean(ZS_CN(i,10:59));
% deltahab(i)=delta_1(i)-delta_2(i);

  
end

idx_resp_hets=idx_hets(find(idx_resp_hets));

figure;imagesc(ZS_CN(idx_resp_hets,:)); colormap('hot');
figure;histogram(delta_1_hets);



std_delta_1_hets=std(delta_1_hets); %%% the standard deviation is simiar to wt... 

idx_delta_pos_hets=find(delta_1_hets>std_delta_1_wt*2);  %%% i am using as a filter 2SD of the wild type delta
idx_delta_neg_hets=find(delta_1_hets<-std_delta_1_wt*2);
idx_delta_hets=union(idx_delta_neg_hets,idx_delta_pos_hets);

histogram(delta_1_hets(idx_delta_hets));

%%% there are no negative ones...

figure;imagesc(ZS_CN(idx_hets(idx_delta_pos_hets),:));  colormap('hot');
figure;plot(mean(ZS_CN(idx_hets(idx_delta_pos_hets),:)));  
  

%%% so i will only use the positive ones

deltahab_2std_hets=delta_1_hets(idx_delta_pos_hets);
Qs = quantile(deltahab_2std_hets,[0.025 0.25 0.50 0.75 0.975]);


idx_delta_pos_hets=idx_hets(idx_delta_pos_hets); %%% to put it in the right indexing format

figure;imagesc(ZS_CN(idx_delta_pos_hets,:));  colormap('hot');
figure;plot(mean(ZS_CN(idx_delta_pos_hets,:)));



%%

figure;
plot(mean(ZS_CN(idx_delta_pos,:)));
hold on;
plot(mean(ZS_CN(idx_delta_pos_fmr1,:)));
hold on;
plot(mean(ZS_CN(idx_delta_pos_hets,:)));
hold off;

%%%% something interesting is that the responses to the 11th loom seem
%%%% stronger in the fmr1 fish!!!

%% per fish


%%% this part is to do it per fish and get the values for the graphs. I
%%% will do it with the hets merged together

fish=unique(idx_Fish);
list5=union(list1,list3);




ROIs_delta_perFish=struct;

%%% this part is for the whole brain. ROIs of each fish that passed the rsq
%%% threshold
for f=1:length(unique(idx_Fish))
    tempfish=find(idx_Fish==fish(f));
    
    
  
    if ismember(fish(f),list2)
        group='fmr1';
        ff=find(ismember(list2,fish(f)));
        temp=intersect(tempfish,idx_delta_pos_fmr1);
    elseif ismember(fish(f),list4)
        group='control';
        ff=find(ismember(list4,fish(f)));
        temp=intersect(tempfish,idx_delta_pos);
    elseif ismember(fish(f),list5)
        group='hets';
        ff=find(ismember(list5,fish(f)));
        temp=intersect(tempfish,idx_delta_pos_hets);
    else 
    end
      
    ROIs_delta_perFish(1).ROIsPerFish.(group)(ff,1)=length(temp);
    ROIs_delta_perFish(2).ROIsPerFish.(group)(ff,1)=length(tempfish);
    
end


%%% now getting the mean and the SD

group = fieldnames(ROIs_delta_perFish(1).ROIsPerFish); %%% mind the order of the fields! 

  
for k=1:3    
ROIs_delta_perFish(1).meanROIsPerFish.mean(k)=mean(ROIs_delta_perFish(1).ROIsPerFish.(group{k}));
ROIs_delta_perFish(1).meanROIsPerFish.SD(k)=std(ROIs_delta_perFish(1).ROIsPerFish.(group{k}));

ROIs_delta_perFish(2).meanROIsPerFish.mean(k)=mean(ROIs_delta_perFish(2).ROIsPerFish.(group{k}));
ROIs_delta_perFish(2).meanROIsPerFish.SD(k)=std(ROIs_delta_perFish(2).ROIsPerFish.(group{k}));

    ROIs_delta_perFish(1).meanROIsPerFish.mean(k)=mean(ROIs_delta_perFish(1).ROIsPerFish.(group{k}));
    ROIs_delta_perFish(1).meanROIsPerFish.SD(k)=std(ROIs_delta_perFish(1).ROIsPerFish.(group{k}));
    
    ROIs_delta_perFish(2).meanROIsPerFish.mean(k)=mean(ROIs_delta_perFish(2).ROIsPerFish.(group{k}));
    ROIs_delta_perFish(2).meanROIsPerFish.SD(k)=std(ROIs_delta_perFish(2).ROIsPerFish.(group{k}));
end


%%% to make the csv files with the means and the SDs. the order is first column is controls, 2nd fmr1
%%% and 3rd column hets
meanROIsDeltaPerFish=[];
meanROIsDeltaPerFish(1,:)=ROIs_CL10_perFish(1).meanROIsPerFish.mean;
meanROIsDeltaPerFish(2,:)=ROIs_CL10_perFish(1).meanROIsPerFish.SD; 


csvwrite(strcat('meanROIsDeltaPerFish','.csv'),meanROIsDeltaPerFish);



%% reclustering

%%%% now this is to cluster each of the fish genotypes that passed the first correlation filter.
%%% I will randomly sample half of the ROIs from the hets to have more or less the same
%%%% amount of ROIs. 

load('s20_fmr1_loomhab_CN_part2.mat','idx_corr');

%%% to also check with different values of the r2
load('rsquare_fmr1loomhab.mat','rsquare_loom');


idx_rsq02=find(rsquare_loom>0.2 & rsquare_loom<1); %%%% it doenst look good at all...

idx_hets_sample= datasample(idx_hets,length(idx_hets)/2,'Replace',false);

Kmeans_perGenotype_postCorr=struct;

options = statset('UseParallel',1);
for i=1:3
  if i==1  
  
      group='fmr1';
      idx_temp=intersect(idx_fmr1,idx_rsq);
    [idxKmeans_ZS_CN_temp Cmap_ZS_CN_temp]=kmeans(ZS_CN(idx_temp,:),100,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    [Model_ZS_CN2,GoodBetas_ZS_CN2]=Test_Regress(Cmap_ZS_CN_temp,Cmap_ZS_CN_temp,idxKmeans_ZS_CN_temp,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

    Kmeans_perGenotype_postCorr.(group).idxKmeans_ZS_CN=idxKmeans_ZS_CN_temp;
    Kmeans_perGenotype_postCorr.(group).Cmap_ZS_CN=Cmap_ZS_CN_temp;
    Kmeans_perGenotype_postCorr.(group).idx=idx_temp;
    
    
  elseif i==2
      
      group='control';
     idx_temp=intersect(idx_wt,idx_rsq);
    [idxKmeans_ZS_CN_temp Cmap_ZS_CN_temp]=kmeans(ZS_CN(idx_temp,:),100,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    [Model_ZS_CN2,GoodBetas_ZS_CN2]=Test_Regress(Cmap_ZS_CN_temp,Cmap_ZS_CN_temp,idxKmeans_ZS_CN_temp,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3
 
     Kmeans_perGenotype_postCorr.(group).idxKmeans_ZS_CN=idxKmeans_ZS_CN_temp;
    Kmeans_perGenotype_postCorr.(group).Cmap_ZS_CN=Cmap_ZS_CN_temp;
    Kmeans_perGenotype_postCorr.(group).idx=idx_temp;
    
  else 
      
      group='hets';
    idx_temp=intersect(idx_hets_sample,idx_rsq);
    [idxKmeans_ZS_CN_temp Cmap_ZS_CN_temp]=kmeans(ZS_CN(idx_temp,:),100,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
    [Model_ZS_CN2,GoodBetas_ZS_CN2]=Test_Regress(Cmap_ZS_CN_temp,Cmap_ZS_CN_temp,idxKmeans_ZS_CN_temp,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

     Kmeans_perGenotype_postCorr.(group).idxKmeans_ZS_CN=idxKmeans_ZS_CN_temp;
    Kmeans_perGenotype_postCorr.(group).Cmap_ZS_CN=Cmap_ZS_CN_temp;
    Kmeans_perGenotype_postCorr.(group).idx=idx_temp;
    
  end  
    
end

%%% here I am sorting the cluster in terms of its representation among the
%%% fish of that group. verygood>9/10, good>3/4 and low<1/2
group=fieldnames(Kmeans_perGenotype_postCorr);
for k=1:3 
    
      idx_temp=Kmeans_perGenotype_postCorr.(group{k}).idx;

for j=1:length(unique(Kmeans_perGenotype_postCorr.(group{k}).idxKmeans_ZS_CN))
    idx_temp2=find(Kmeans_perGenotype_postCorr.(group{k}).idxKmeans_ZS_CN==j);
    fish_temp=idx_Fish(idx_temp(idx_temp2));

    Kmeans_perGenotype_postCorr.(group{k}).fishN.clusters(j)=length(unique(fish_temp));
    end
    verygood=find(Kmeans_perGenotype_postCorr.(group{k}).fishN.clusters>9*(length(unique(fish_temp))/10));
    good=find(Kmeans_perGenotype_postCorr.(group{k}).fishN.clusters>3*(length(unique(fish_temp))/4));
    low=find(Kmeans_perGenotype_postCorr.(group{k}).fishN.clusters<length(unique(fish_temp))/2);
    Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_verygood=Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN(verygood,:);
    Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_good=Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN(good,:);
    Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_low=Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN(low,:);
    

end



%% testing their similarity with a correlation

%%% not sure if this will really be useful. 

%%% i am also testing a hierarchical clustering based on the correlation.
%%% is not looking  bad... but I also separated them by the clusters that are represented in most fish
%%% the results are nice but not very different that my original Kmeans...
%%% at least when analyzed the 3 genotypes together. 


start=1;
for k=1:3
 
    temp_corr=corrcoef(Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_verygood');
    Kmeans_perGenotype_postCorr.(group{k}).corr_verygood=temp_corr;
    position=[start:start+length(temp_corr)-1];
    Kmeans_perGenotype_postCorr.all.Cmap_ZS_CN_verygood(position,:)=Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_verygood;
        
      start=max(position)+1;
      
end    




temp_corr=corrcoef(Kmeans_perGenotype_postCorr.all.Cmap_ZS_CN_verygood');
Kmeans_perGenotype_postCorr.all.corr_verygood=temp_corr;

Z = linkage(temp_corr,'complete','correlation');

dendrogram(Z);

T = cluster(Z,'cutoff',0.2,'criterion','distance');
unique(T)

dendrogram(Z,0,'colorthreshold',0.2);

figure; counter=1;
for i=1:length(unique(T))

    subplot(4,4,counter);
    
    idx_temp=find(T==i);
    if length(idx_temp)>1
    plot(mean(Kmeans_perGenotype_postCorr.all.Cmap_ZS_CN_verygood(idx_temp,:)));
    else
     plot(Kmeans_perGenotype_postCorr.all.Cmap_ZS_CN_verygood(idx_temp,:));   
    end
    
    counter=counter+1;
end


%%

%%%% I will do the same but for each genotype. 
%%% I am not sure that what I am finding is interesting... 


for k=1:3
temp_corr=corrcoef(Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_verygood');
Kmeans_perGenotype_postCorr.(group{k}).corr_verygood=temp_corr;

Z = linkage(temp_corr,'complete','correlation');

T = cluster(Z,'cutoff',0.25,'criterion','distance')
unique(T)

figure;dendrogram(Z,0,'colorthreshold',0.25)

figure; counter=1;
for i=1:length(unique(T))

    subplot(4,4,counter);
    
    idx_temp=find(T==i);
    if length(idx_temp)>1
    plot(mean(Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_verygood(idx_temp,:)));
    else
     plot(Kmeans_perGenotype_postCorr.(group{k}).Cmap_ZS_CN_verygood(idx_temp,:));   
    end
    
    counter=counter+1;
end
end



