%%%%% this is part3 of the fmr1loomhab analysis

%%% I am mostly getting some things to make some figures for the M3
%%% I am planning to get the means, the max responses, the distribution of
%%% the max responses (histograms) and I will try to get something with the
%%% localization of ROIs.

%%%%% NOTE: I need to confirm with lena but I think list 2 is the fmr1 fish
%%%%% list 4 is the wildtypes and the others are the hets


load('s20_fmr1_loomhab_CN.mat','MatFiles','ZS_CN');


load('s20_fmr1_loomhab_CN_post_rsq01.mat','idx_rsq1');

load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');


%%% to get the clusters from the Kmeans done to ALL the ROIs
load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN');

load('s20_good_idx_Fish.mat','idx_Fish');
idx_Fish_cat=categorical(idx_Fish);

%%

load('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','idx_rsq');

load('s20_fmr1_loomhab_CN_part2_High_corr_Nb.mat','High_corr_Nb');

load('fmr1loomhab_BrainRegNclean.mat','PerBrainRegions','RegionList','ROI_temp2','idx_rsq_cleaned');

%%% cleaning the clasification of the clusters. 
idx_clean=ismember(idx_rsq,idx_rsq_cleaned);
idx_clean=find(idx_clean);

High_corr_Nb=High_corr_Nb(idx_clean);

%%% check if it worked
figure;counter=1;
for i=1:length(unique(High_corr_Nb))
    
    idx_temp=idx_rsq_cleaned(find(High_corr_Nb==i));
    
    subplot(3,4,counter)
    plot(mean(ZS_CN(idx_temp,:)))
    counter=counter+1;
end

%%% it worked

 %%%% i also need to generate the list of fish saved in the
 %%%% fmr1loomhab_lists.m file

 %%
 %%% this is to make a variable with what I need for my analysis in a
%%% shorter version. 
 
GoodClust_goodmembers_full=[];
for i=1:length(unique(High_corr_Nb))    
    idx_temp=idx_rsq_cleaned(find(High_corr_Nb==i));
    
    GoodClust_goodmembers_full(i).ZS=ZS_CN(idx_temp,:);
    GoodClust_goodmembers_full(i).idx=idx_temp;
    GoodClust_goodmembers_full(i).mean=mean(GoodClust_goodmembers_full(i).ZS,1);
    GoodClust_goodmembers_full(i).STD=std(GoodClust_goodmembers_full(i).ZS,1,1);       
end


%%
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 600]);
rows=length(unique(High_corr_Nb));counter=1;
for i=1:length(unique(High_corr_Nb))
    
    subplot(rows,4,counter);plot(GoodClust_goodmembers_full(i).mean); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(GoodClust_goodmembers_full(i).ZS,[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(GoodClust_goodmembers_full(i).idx)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_cat(GoodClust_goodmembers_full(i).idx));% ax=gca; ax.FontSize=4;%%% for the fish location
    counter=counter+4;
    
end


%%
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 600]);
rows=length(unique(High_corr_Nb))-2;counter=1;
for i=3:length(unique(High_corr_Nb))
    
    subplot(rows,2,counter,'replace');plot(GoodClust_goodmembers_full(i).mean); ylim([-2 10]) %%%to plot the mean
    subplot(rows,2,counter+1,'replace');imagesc(GoodClust_goodmembers_full(i).ZS,[-2 10]); colormap('hot')%%%for the raster plot
    %subplot(rows,4,counter+2);histogram(idx_Plane(GoodClust_goodmembers_full(i).idx)); %%%for the plane location   
    %subplot(rows,3,counter+2,'replace');histogram(idx_Fish_cat(GoodClust_goodmembers_full(i).idx));% ax=gca; ax.FontSize=4;%%% for the fish location
    counter=counter+2;
    
end

%%% this part is just to check at the fast hab ones. it seems that I only
%%% have a sharp and broads. but this lat got separated...

figure;

for i=[3 7 8];
plot(GoodClust_goodmembers_full(i).mean(1,60:100)); ylim([-2 10]);
hold on;
end



%%
%%% to check the means


idx_temp1=ismember(idx_Fish,list1);
idx_temp1=find(idx_temp1);
idx_temp2=ismember(idx_Fish,list2);
idx_temp2=find(idx_temp2);
idx_temp3=ismember(idx_Fish,list3);
idx_temp3=find(idx_temp3);
idx_temp4=ismember(idx_Fish,list4);
idx_temp4=find(idx_temp4);


%%% this is to merge the 2 het groups
idx_temp5=union(idx_temp1,idx_temp3);


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 600]);
counter=1;
for i=3:length(unique(High_corr_Nb)) %%%% i am avoiding the first and 2nd clusters. the first is not representative and the 2nd one looks like artifact. 
     subplot(8,1,counter);
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp5);temp=find(temp); %%% this is for the hets together
   plot(mean(GoodClust_goodmembers_full(i).ZS(temp,:)));ylim([-2 10]) %%%to plot the mean
   hold on;  
     
     
%    temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp1);temp=find(temp);
%    plot(mean(GoodClust_goodmembers_full(i).ZS(temp,:)));ylim([-2 10]) %%%to plot the mean
%    hold on;
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp2);temp=find(temp);
   plot(mean(GoodClust_goodmembers_full(i).ZS(temp,:))); ylim([-2 10])%%%to plot the mean
   hold on;
   
%    temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp3);temp=find(temp);
%    plot(mean(GoodClust_goodmembers_full(i).ZS(temp,:)));ylim([-2 10]) %%%to plot the mean
%    hold on;
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp4);temp=find(temp);
   plot(mean(GoodClust_goodmembers_full(i).ZS(temp,:)));ylim([-2 10]) %%%to plot the mean
   hold on;

counter=counter+1;
end

%%% to get the means in cvs files
%%% in the order of the figures [3 7 8 6 5 4 9 10]

fmr1_loomhab_CL10_means=[];
counter=1;
for i=[3 7 8 6 5 4 9 10]  %%% i am not taking the first 2 clusters
      
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp2);temp=find(temp);
   fmr1_loomhab_CL10_means(counter,:)=mean(GoodClust_goodmembers_full(i).ZS(temp,:)); %%%to plot the mean
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp4);temp=find(temp);
   fmr1_loomhab_CL10_means(counter+1,:)=mean(GoodClust_goodmembers_full(i).ZS(temp,:)); %%%to plot the mean
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp5);temp=find(temp); %%% this is for the hets together
   fmr1_loomhab_CL10_means(counter+2,:)=mean(GoodClust_goodmembers_full(i).ZS(temp,:)); %%%to plot the mean
   
    counter=counter+6;
end

 csvwrite(strcat('fmr1_loomhab_CL10_means','_Allclust','.csv'),fmr1_loomhab_CL10_means);
  





%%

%%% this part is to check at the max responses, but per cluster. I will not
%%% use cluster 1 as is mostly 2 fish only, neither the inhibited responses

%%%% this is based in a part where gilles tried to look at the distribution
%%%% of the max responses of the 3 first looms. he was filtering by neurons
%%%% that responsed 3 "SD" above but i am not doing it anymore cause i
%%%% filtered with single linear regressions


for i=1:length(unique(High_corr_Nb)) 
    
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp1);temp=find(temp);

idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(temp,60:80),[],2);
GoodClust_goodmembers_full(i).Responses_strength{1,1}=idx_temp_res;  
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,110:130),[],2);
GoodClust_goodmembers_full(i).Responses_strength{1,2}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,145:165),[],2);
GoodClust_goodmembers_full(i).Responses_strength{1,3}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,185:205),[],2);%Loom 4
GoodClust_goodmembers_full(i).Responses_strength{1,4}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,220:240),[],2);%Loom 5
GoodClust_goodmembers_full(i).Responses_strength{1,5}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,500:550),[],2);%Loom 11
GoodClust_goodmembers_full(i).Responses_strength{1,6}=idx_temp_res2;
   
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp2);temp=find(temp);
   
idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(temp,60:80),[],2);
GoodClust_goodmembers_full(i).Responses_strength{2,1}=idx_temp_res;  
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,110:130),[],2);
GoodClust_goodmembers_full(i).Responses_strength{2,2}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,145:165),[],2);
GoodClust_goodmembers_full(i).Responses_strength{2,3}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,185:205),[],2);%Loom 4
GoodClust_goodmembers_full(i).Responses_strength{2,4}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,220:240),[],2);%Loom 5
GoodClust_goodmembers_full(i).Responses_strength{2,5}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,500:550),[],2);%Loom 11
GoodClust_goodmembers_full(i).Responses_strength{2,6}=idx_temp_res2;
   
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp3);temp=find(temp);
   
idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(temp,60:80),[],2);
GoodClust_goodmembers_full(i).Responses_strength{3,1}=idx_temp_res;  
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,110:130),[],2);
GoodClust_goodmembers_full(i).Responses_strength{3,2}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,145:165),[],2);
GoodClust_goodmembers_full(i).Responses_strength{3,3}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,185:205),[],2);%Loom 4
GoodClust_goodmembers_full(i).Responses_strength{3,4}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,220:240),[],2);%Loom 5
GoodClust_goodmembers_full(i).Responses_strength{3,5}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,500:550),[],2);%Loom 11
GoodClust_goodmembers_full(i).Responses_strength{3,6}=idx_temp_res2;
   
   
temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp4);temp=find(temp);

idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(temp,60:80),[],2);
GoodClust_goodmembers_full(i).Responses_strength{4,1}=idx_temp_res;  
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,110:130),[],2);
GoodClust_goodmembers_full(i).Responses_strength{4,2}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,145:165),[],2);
GoodClust_goodmembers_full(i).Responses_strength{4,3}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,185:205),[],2);%Loom 4
GoodClust_goodmembers_full(i).Responses_strength{4,4}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,220:240),[],2);%Loom 5
GoodClust_goodmembers_full(i).Responses_strength{4,5}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,500:550),[],2);%Loom 11
GoodClust_goodmembers_full(i).Responses_strength{4,6}=idx_temp_res2;
  

end

%%%% this is for the hets together
for i=1:length(unique(High_corr_Nb)) 
    
   temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp5);temp=find(temp);

idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(temp,60:80),[],2);
GoodClust_goodmembers_full(i).Responses_strength{5,1}=idx_temp_res;  
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,110:130),[],2);
GoodClust_goodmembers_full(i).Responses_strength{5,2}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,145:165),[],2);
GoodClust_goodmembers_full(i).Responses_strength{5,3}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,185:205),[],2);%Loom 4
GoodClust_goodmembers_full(i).Responses_strength{5,4}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,220:240),[],2);%Loom 5
GoodClust_goodmembers_full(i).Responses_strength{5,5}=idx_temp_res2;
idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,500:550),[],2);%Loom 11
GoodClust_goodmembers_full(i).Responses_strength{5,6}=idx_temp_res2;

end



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 900]);
edges=[0:0.25:15];
rows=length(unique(High_corr_Nb))-2;counter=1;
for i=2:length(unique(High_corr_Nb))-1 %%%% i am avoiding the first and last clusters. the first is not representative and the last is the inhibited one. 
 
subplot(rows,6,counter);
for j=[5 2 4] %1:4
   histogram(GoodClust_goodmembers_full(i).Responses_strength{j,1},edges,'normalization','probability');
    hold on; 
end

subplot(rows,6,counter+1);

for j=[5 2 4] %1:4
   histogram(GoodClust_goodmembers_full(i).Responses_strength{j,2},edges,'normalization','probability');
    hold on; 
end

subplot(rows,6,counter+2);

for j=[5 2 4] %1:4
   histogram(GoodClust_goodmembers_full(i).Responses_strength{j,3},edges,'normalization','probability');
    hold on; 
end

subplot(rows,6,counter+3);

for j=[5 2 4] %1:4
   histogram(GoodClust_goodmembers_full(i).Responses_strength{j,4},edges,'normalization','probability');
    hold on; 
end

subplot(rows,6,counter+4);

for j=[5 2 4] %1:4
   histogram(GoodClust_goodmembers_full(i).Responses_strength{j,5},edges,'normalization','probability');
    hold on; 
end

subplot(rows,6,counter+5);

for j=[5 2 4] %1:4
   histogram(GoodClust_goodmembers_full(i).Responses_strength{j,6},edges,'normalization','probability');
    hold on; 
end


counter=counter+6;

end






%% per fish
%%% this part is to do it per fish and get the values for the graphs. I
%%% will do it with the hets merged together

fish=unique(idx_Fish);
list5=union(list1,list3);
edges=[0:0.25:15];
%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know why.


for f=1:length(unique(idx_Fish))
    tempfish=find(idx_Fish==fish(f));
    
  for i=1:length(unique(High_corr_Nb))
      
      temp=ismember(GoodClust_goodmembers_full(i).idx,tempfish);temp=find(temp);
      
    if ismember(fish(f),list2)
        group='fmr1';
        ff=find(ismember(list2,fish(f)));
    elseif ismember(fish(f),list4)
        group='control';
        ff=find(ismember(list4,fish(f)));
    elseif ismember(fish(f),list5)
        group='hets';
        ff=find(ismember(list5,fish(f)));
    else 
    end
     
    
    idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(temp,60:80),[],2);
    GoodClust_goodmembers_full(i).MaxPerFish.(group){ff,1}=idx_temp_res;
    [GoodClust_goodmembers_full(i).MaxCountPerFish.(group).loom1.count(ff,:),~]=histcounts(idx_temp_res,edges);
    
    
    idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,110:130),[],2);
    GoodClust_goodmembers_full(i).MaxPerFish.(group){ff,2}=idx_temp_res2;
    [GoodClust_goodmembers_full(i).MaxCountPerFish.(group).loom2.count(ff,:),~]=histcounts(idx_temp_res2,edges);
     
    idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,145:165),[],2);
    GoodClust_goodmembers_full(i).MaxPerFish.(group){ff,3}=idx_temp_res2;
    [GoodClust_goodmembers_full(i).MaxCountPerFish.(group).loom3.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
    idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,185:205),[],2);%Loom 4
    GoodClust_goodmembers_full(i).MaxPerFish.(group){ff,4}=idx_temp_res2;
    [GoodClust_goodmembers_full(i).MaxCountPerFish.(group).loom4.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
    idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,220:240),[],2);%Loom 5
    GoodClust_goodmembers_full(i).MaxPerFish.(group){ff,5}=idx_temp_res2;
    [GoodClust_goodmembers_full(i).MaxCountPerFish.(group).loom5.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
    idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(temp,500:550),[],2);%Loom 11
    GoodClust_goodmembers_full(i).MaxPerFish.(group){ff,6}=idx_temp_res2; 
    [GoodClust_goodmembers_full(i).MaxCountPerFish.(group).loom6.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
end
end

for f=1:length(unique(idx_Fish))
for i=1:length(unique(High_corr_Nb))

 if ismember(fish(f),list2)
        group='fmr1';
        ff=find(ismember(list2,fish(f)));
    elseif ismember(fish(f),list4)
        group='control';
        ff=find(ismember(list4,fish(f)));
    elseif ismember(fish(f),list5)
        group='hets';
        ff=find(ismember(list5,fish(f)));
    else 
 end   
  
for loom=1:6    
GoodClust_goodmembers_full(i).MaxCountPerFish.(group).(strcat('loom',num2str(loom))).mean(1,:)=mean(GoodClust_goodmembers_full(i).MaxCountPerFish.(group).(strcat('loom',num2str(loom))).count);
end
end
end

save('s20_fmr1_loomhab_CN_part3.mat','GoodClust_goodmembers_full','idx_temp1','idx_temp2','idx_temp3','idx_temp4','idx_temp5');

%%% this part is to get some of usual clusters and get the means of the
%%% distribution to plot them. 

group = fieldnames(GoodClust_goodmembers_full(i).MaxCountPerFish);
counter=1;

 figure;
for i=2:length(unique(High_corr_Nb))-1
   BinnedRespStrLooms_rsq_perFish=[];
   
for loom=1:6
subplot(8,6,counter);
    for k=1:3
     
     plot(GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean);
     hold on;

BinnedRespStrLooms_rsq_perFish(k,:)=GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean;
  
    end
     hold off;
     counter=counter+1;
   
%csvwrite(strcat('BinnedRespStrLooms_rsq_perFish_','clust',num2str(i),'.csv'),BinnedRespStrLooms_rsq_perFish');
 
end
end

%%% this part is to get some of usual clusters and get the means of the
%%% distribution to plot them in the right order. 

group = fieldnames(GoodClust_goodmembers_full(i).MaxCountPerFish);
counter=1;

 figure;
for i=[3 7 8 6 5 4 9 10] %%% I am not taking the first 2 clusters. 
   
for loom=1:6
subplot(8,6,counter);
    for k=1:3
     
     plot(GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean); 
     hold on;
  
    end
     hold off;
     counter=counter+1;
    
end
end


%%

%%% i still need to work in this part
 
%%%% now the max responses of the first 3 looms... in box plots and changed
%%%% histogram

idx_temp_res=max(ZS_CN(idx_rsq_cleaned,60:80),[],2);
Responses_strength{1}=idx_temp_res(idx_temp_res>threshold);
idx_temp_res2=max(ZS_CN(idx_rsq_cleaned,110:130),[],2);
Responses_strength{2}=idx_temp_res2(idx_temp_res>threshold);
idx_temp_res2=max(ZS_CN(idx_rsq_cleaned,145:165),[],2);
Responses_strength{3}=idx_temp_res2(idx_temp_res>threshold);

idx_Fish_rsq=idx_Fish(idx_rsq_cleaned);
idx_Fish_rsq=idx_Fish_rsq(idx_temp_res>threshold);
idx_temp1=ismember(idx_Fish_rsq,list1);
idx_temp1=find(idx_temp1);
idx_temp2=ismember(idx_Fish_rsq,list2);
idx_temp2=find(idx_temp2);
idx_temp3=ismember(idx_Fish_rsq,list3);
idx_temp3=find(idx_temp3);
idx_temp4=ismember(idx_Fish_rsq,list4);
idx_temp4=find(idx_temp4);


Responses_strength_forboxplot={};
Responses_strength_forboxplot{1}=[Responses_strength{1}(idx_temp1)', Responses_strength{1}(idx_temp2)',Responses_strength{1}(idx_temp3)',Responses_strength{1}(idx_temp4)'];
Responses_strength_forboxplot{2}=[Responses_strength{2}(idx_temp1)', Responses_strength{2}(idx_temp2)',Responses_strength{2}(idx_temp3)',Responses_strength{2}(idx_temp4)'];
Responses_strength_forboxplot{3}=[Responses_strength{3}(idx_temp1)', Responses_strength{3}(idx_temp2)',Responses_strength{3}(idx_temp3)',Responses_strength{3}(idx_temp4)'];

grp = [ones(1,length(idx_temp1)),2*ones(1,length(idx_temp2)),3*ones(1,length(idx_temp3)),4*ones(1,length(idx_temp4))];


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 900]);

for i=1:3
subplot(1,3,i);
boxplot(Responses_strength_forboxplot{i},grp); ylim([-2 20]);

end


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 900]);
edges=[0:0.5:15];
subplot(1,3,1);
histogram(Responses_strength{1}(idx_temp1),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{1}(idx_temp2),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{1}(idx_temp3),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{1}(idx_temp4),edges,'normalization','probability','DisplayStyle', 'stairs');
subplot(1,3,2);
histogram(Responses_strength{2}(idx_temp1),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{2}(idx_temp2),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{2}(idx_temp3),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{2}(idx_temp4),edges,'normalization','probability','DisplayStyle', 'stairs');
subplot(1,3,3);
histogram(Responses_strength{3}(idx_temp1),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{3}(idx_temp2),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{3}(idx_temp3),edges,'normalization','probability','DisplayStyle', 'stairs');hold on;histogram(Responses_strength{3}(idx_temp4),edges,'normalization','probability','DisplayStyle', 'stairs');


%% ROIs whole brain count

%%%% this part is to do a rough count of the ROIs in the whole brain per
%%%% cluster.  I need to normilize it by fish. 

ROIs_CL10_perFish=struct;

%%% this part is for the whole brain. ROIs of each fish that passed the rsq
%%% threshold
for f=1:length(unique(idx_Fish))
    tempfish=find(idx_Fish==fish(f));
    
    temp=intersect(tempfish,idx_rsq_cleaned);
  
    if ismember(fish(f),list2)
        group='fmr1';
        ff=find(ismember(list2,fish(f)));
    elseif ismember(fish(f),list4)
        group='control';
        ff=find(ismember(list4,fish(f)));
    elseif ismember(fish(f),list5)
        group='hets';
        ff=find(ismember(list5,fish(f)));
    else 
    end
      
    ROIs_CL10_perFish(11).ROIsPerFish.(group)(ff,1)=length(temp);
    ROIs_CL10_perFish(12).ROIsPerFish.(group)(ff,1)=length(tempfish);
end

%%% this is to do it per clusters
for f=1:length(unique(idx_Fish))
    tempfish=find(idx_Fish==fish(f));
    
  for i=1:length(unique(High_corr_Nb))
      
      temp=ismember(GoodClust_goodmembers_full(i).idx,tempfish);temp=find(temp);
      
    if ismember(fish(f),list2)
        group='fmr1';
        ff=find(ismember(list2,fish(f)));
    elseif ismember(fish(f),list4)
        group='control';
        ff=find(ismember(list4,fish(f)));
    elseif ismember(fish(f),list5)
        group='hets';
        ff=find(ismember(list5,fish(f)));
    else 
    end
     
    
    ROIs_CL10_perFish(i).ROIsPerFish.(group)(ff,1)=length(temp);
    
    
end
end


%%% now getting the mean and the SD
group = fieldnames(ROIs_CL10_perFish(i).ROIsPerFish); %%% mind the order of the fields! 

   
for i=1:length(unique(High_corr_Nb))

  
for k=1:3    
ROIs_CL10_perFish(i).meanROIsPerFish.mean(k)=mean(ROIs_CL10_perFish(i).ROIsPerFish.(group{k}));
ROIs_CL10_perFish(i).meanROIsPerFish.SD(k)=std(ROIs_CL10_perFish(i).ROIsPerFish.(group{k}));


    ROIs_CL10_perFish(11).meanROIsPerFish.mean(k)=mean(ROIs_CL10_perFish(11).ROIsPerFish.(group{k}));
    ROIs_CL10_perFish(11).meanROIsPerFish.SD(k)=std(ROIs_CL10_perFish(11).ROIsPerFish.(group{k}));
    
    ROIs_CL10_perFish(12).meanROIsPerFish.mean(k)=mean(ROIs_CL10_perFish(12).ROIsPerFish.(group{k}));
    ROIs_CL10_perFish(12).meanROIsPerFish.SD(k)=std(ROIs_CL10_perFish(12).ROIsPerFish.(group{k}));
end
end


%%% to make the csv files with the means and the SDs. the order is first column is controls, 2nd fmr1
%%% and 3rd column hets

meanROIsPerFish=[];
counter=1;
for i=2:length(unique(High_corr_Nb)) %%% I am not taking the first cluster
   
        
meanROIsPerFish(counter,:)=ROIs_CL10_perFish(i).meanROIsPerFish.mean;
meanROIsPerFish(counter+1,:)=ROIs_CL10_perFish(i).meanROIsPerFish.SD;  


counter=counter+4; %%% this is to make some espace between cluster data

 
end
%%% and for the whole brain rsq ROIs
meanROIsPerFish(counter,:)=ROIs_CL10_perFish(11).meanROIsPerFish.mean;
meanROIsPerFish(counter+1,:)=ROIs_CL10_perFish(11).meanROIsPerFish.SD;  

counter=counter+4; %%% this is to make some espace between cluster data
%%% and for the whole brain ROIs including non respondent

meanROIsPerFish(counter,:)=ROIs_CL10_perFish(12).meanROIsPerFish.mean;
meanROIsPerFish(counter+1,:)=ROIs_CL10_perFish(12).meanROIsPerFish.SD; 


csvwrite(strcat('meanROIsPerFish','_Allclust','.csv'),meanROIsPerFish);




