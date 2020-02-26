
%%%%% this script was an extra of part3analysis_fmr1loomhab.m just to look
%%%%% at the  loom responses to all the rsq ROIs and not per cluster. I
%%%%% dont think is that useful. 

MaxPerFish_rsq=[];


fish=unique(idx_Fish);
list5=union(list1,list3);
edges=[0:0.25:15];
%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know why.
for f=1:length(unique(idx_Fish))
    tempfish=find(idx_Fish==fish(f));
      
      temp=intersect(idx_rsq,tempfish);
      
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
     
    
    idx_temp_res=max(ZS_CN(temp,60:80),[],2);
    MaxPerFish_rsq.MaxPerFish.(group){ff,1}=idx_temp_res;
    [MaxPerFish_rsq.MaxCountPerFish.(group).loom1.count(ff,:),~]=histcounts(idx_temp_res,edges);
    
    
    idx_temp_res2=max(ZS_CN(temp,110:130),[],2);
    MaxPerFish_rsq.MaxPerFish.(group){ff,2}=idx_temp_res2;
    [MaxPerFish_rsq.MaxCountPerFish.(group).loom2.count(ff,:),~]=histcounts(idx_temp_res2,edges);
     
    idx_temp_res2=max(ZS_CN(temp,145:165),[],2);
    MaxPerFish_rsq.MaxPerFish.(group){ff,3}=idx_temp_res2;
    [MaxPerFish_rsq.MaxCountPerFish.(group).loom3.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
    idx_temp_res2=max(ZS_CN(temp,185:205),[],2);%Loom 4
    MaxPerFish_rsq.MaxPerFish.(group){ff,4}=idx_temp_res2;
    [MaxPerFish_rsq.MaxCountPerFish.(group).loom4.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
    idx_temp_res2=max(ZS_CN(temp,220:240),[],2);%Loom 5
    MaxPerFish_rsq.MaxPerFish.(group){ff,5}=idx_temp_res2;
    [MaxPerFish_rsq.MaxCountPerFish.(group).loom5.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
    idx_temp_res2=max(ZS_CN(temp,500:550),[],2);%Loom 11
    MaxPerFish_rsq.MaxPerFish.(group){ff,6}=idx_temp_res2; 
    [MaxPerFish_rsq.MaxCountPerFish.(group).loom6.count(ff,:),~]=histcounts(idx_temp_res2,edges);
    
end

for f=1:length(unique(idx_Fish))

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
MaxPerFish_rsq.MaxCountPerFish.(group).(strcat('loom',num2str(loom))).mean(1,:)=mean(MaxPerFish_rsq.MaxCountPerFish.(group).(strcat('loom',num2str(loom))).count);
end
end


%%% this part is to get some of usual clusters and get the means of the
%%% distribution to plot them. 

group = fieldnames(MaxPerFish_rsq.MaxCountPerFish);
counter=1;

 figure;

   Whole_BinnedRespStrLooms_rsq_perFish=[];
   
for loom=1:6
subplot(1,6,counter);
    for k=1:3
     
     plot(MaxPerFish_rsq.MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean);
     hold on;

Whole_BinnedRespStrLooms_rsq_perFish(k,:)=MaxPerFish_rsq.MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean;
  
    end
     hold off;
     counter=counter+1;
   
%csvwrite(strcat('Whole_BinnedRespStrLooms_rsq_perFish','.csv'),Whole_BinnedRespStrLooms_rsq_perFish');
 
end


%%
edges=[0.25:0.25:15];
group = fieldnames(MaxPerFish_rsq.MaxCountPerFish);
counter=1;

 figure;

   Whole_BinnedRespStrLooms_rsq_perFish=[];
   
for loom=[1 2 3 6];
subplot(1,4,counter);
    for k=1:3
     
     %plot(MaxPerFish_rsq.MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean);
     

     temp_area=trapz(edges,MaxPerFish_rsq.MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean);
     temp_mean=MaxPerFish_rsq.MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean;
     temp_mean_norm=temp_mean./temp_area;
     %trapz(edges, temp_mean_norm) %%% to check if it work. the answer
     %should be 1.
     if k==1
       alpha=15;  
     elseif k==2
         alpha=0;
     else
         alpha=10;
     end
     
     plot(edges,temp_mean_norm,'linewidth',1.5,'Color',[0 0 0]+0.05*alpha);
     hold on;
     
     
     
     
Whole_BinnedRespStrLooms_rsq_perFish(k,:)=MaxPerFish_rsq.MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean;
  
    end
     hold off;
     counter=counter+1;
   
%csvwrite(strcat('Whole_BinnedRespStrLooms_rsq_perFish','.csv'),Whole_BinnedRespStrLooms_rsq_perFish');
 
end