

%%

%%% this part is to check at the max responses, but per cluster and all the genotypes together 
%%% I am looking at the >2SD higher end of the ROIs  
%%%  and >1SD. this doesnt look that interesting.
%%%  and >3SD. this one is interesting but not as much as the >2SD.

for i=1:length(unique(High_corr_Nb)) 
     
   %%% for fmr1
   %temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp2);temp=find(temp); 
   
        
%%% for controls  
%temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp4);temp=find(temp);  


%%%% this is for the hets together  
%temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp5);temp=find(temp);
   
idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(:,60:80),[],2);
temp_std=std(idx_temp_res);
temp_mean=mean(idx_temp_res);
idx_2SD_resp=find(idx_temp_res>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_2SD{1,1}=GoodClust_goodmembers_full(i).idx(idx_2SD_resp);
GoodClust_goodmembers_full(i).Responses_strength{6,1}=idx_temp_res;


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,110:130),[],2);
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_2SD_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_2SD{1,2}=GoodClust_goodmembers_full(i).idx(idx_2SD_resp);
GoodClust_goodmembers_full(i).Responses_strength{6,2}=idx_temp_res2;


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,145:165),[],2);
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_2SD_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_2SD{1,3}=GoodClust_goodmembers_full(i).idx(idx_2SD_resp);
GoodClust_goodmembers_full(i).Responses_strength{6,3}=idx_temp_res2;


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,185:205),[],2);%Loom 4
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_2SD_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_2SD{1,4}=GoodClust_goodmembers_full(i).idx(idx_2SD_resp);
GoodClust_goodmembers_full(i).Responses_strength{6,4}=idx_temp_res2;


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,220:240),[],2);%Loom 5
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_2SD_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_2SD{1,5}=GoodClust_goodmembers_full(i).idx(idx_2SD_resp);
GoodClust_goodmembers_full(i).Responses_strength{6,5}=idx_temp_res2;


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,500:550),[],2);%Loom 11
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_2SD_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_2SD{1,6}=GoodClust_goodmembers_full(i).idx(idx_2SD_resp);
GoodClust_goodmembers_full(i).Responses_strength{6,6}=idx_temp_res2;


end

%%% to plot them and see how many ROIs per genotype

for i=[3 7 8 6]
    
    temp2=intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp2);
    temp5=intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp5);
    temp4=intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp4);
    
figure;
scatter(ROI_temp2(GoodClust_goodmembers_full(i).idx,2),ROI_temp2(GoodClust_goodmembers_full(i).idx,1),10,'filled');
hold on;
scatter(ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp2),2),ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp2),1),10,'filled');
hold on;
scatter(ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp5),2),ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp5),1),10,'filled');
hold on;
scatter(ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp4),2),ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_2SD{1,2},idx_temp4),1),10,'filled');

title(strcat('clust_',num2str(i),'_fmr1=',num2str(length(temp2)/11),'_hets=',num2str(length(temp5)/20),'_wt=',num2str(length(temp4)/10)));
end

%%

%%% this part is to check at the max responses, but per cluster and all the genotypes together 
%%% I am looking at and arbitrary threshold for the 2nd loom based on the graphs

%%% it gave me similar results than the 2SD...

for i=1:length(unique(High_corr_Nb)) 
     
   %%% for fmr1
   %temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp2);temp=find(temp); 
   
        
%%% for controls  
%temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp4);temp=find(temp);  


%%%% this is for the hets together  
%temp=ismember(GoodClust_goodmembers_full(i).idx,idx_temp5);temp=find(temp);

if i==1
    thresh=2.25;
elseif i==2
    thresh=2;
elseif i==3
    thresh=3;
elseif i==4
    thresh=3.5;
elseif i==5
    thresh=3.75;
else     
end
   
idx_temp_res=max(GoodClust_goodmembers_full(i).ZS(:,60:80),[],2);
temp_std=std(idx_temp_res);
temp_mean=mean(idx_temp_res);
idx_thresh_resp=find(idx_temp_res>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_thresh{1,1}=GoodClust_goodmembers_full(i).idx(idx_thresh_resp);


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,110:130),[],2);
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_thresh_resp=find(idx_temp_res2>thresh);
GoodClust_goodmembers_full(i).idx_thresh{1,2}=GoodClust_goodmembers_full(i).idx(idx_thresh_resp);


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,145:165),[],2);
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_thresh_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_thresh{1,3}=GoodClust_goodmembers_full(i).idx(idx_thresh_resp);


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,185:205),[],2);%Loom 4
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_thresh_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_thresh{1,4}=GoodClust_goodmembers_full(i).idx(idx_thresh_resp);


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,220:240),[],2);%Loom 5
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_thresh_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_thresh{1,5}=GoodClust_goodmembers_full(i).idx(idx_thresh_resp);


idx_temp_res2=max(GoodClust_goodmembers_full(i).ZS(:,500:550),[],2);%Loom 11
temp_std=std(idx_temp_res2);
temp_mean=mean(idx_temp_res2);
idx_thresh_resp=find(idx_temp_res2>temp_mean+2*temp_std);
GoodClust_goodmembers_full(i).idx_thresh{1,6}=GoodClust_goodmembers_full(i).idx(idx_thresh_resp);


end

%%% to plot them and see how many ROIs per genotype

for i=[3 7 8 6]
    
    temp2=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp2);
    temp5=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp5);
    temp4=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp4);
    
figure;
scatter(ROI_temp2(GoodClust_goodmembers_full(i).idx,2),ROI_temp2(GoodClust_goodmembers_full(i).idx,1),10,'filled');
hold on;
scatter(ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp2),2),ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp2),1),10,'filled');
hold on;
scatter(ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp5),2),ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp5),1),10,'filled');
hold on;
scatter(ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp4),2),ROI_temp2(intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp4),1),10,'filled');

title(strcat('clust_',num2str(i),'_fmr1=',num2str(length(temp2)/11),'_hets=',num2str(length(temp5)/20),'_wt=',num2str(length(temp4)/10)));
end


%%% with the zbrain mask
%%% 2nd loom
for i=[3 7 8 6]
    
    temp2=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp2);
    temp5=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp5);
    temp4=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,2},idx_temp4);
    
figure; patch(brain3D,'EdgeColor','none','FaceAlpha',0.05);hold on;
% scatter3(ROI_temp2(GoodClust_goodmembers_full(i).idx,1),ROI_temp2(GoodClust_goodmembers_full(i).idx,2),ROI_temp2(GoodClust_goodmembers_full(i).idx,3),10,'filled');
% hold on;
scatter3(ROI_temp2(temp2,1),ROI_temp2(temp2,2),ROI_temp2(temp2,3),10,'filled');
hold on;
scatter3(ROI_temp2(temp5,1),ROI_temp2(temp5,2),ROI_temp2(temp5,3),10,'filled');
hold on;
scatter3(ROI_temp2(temp4,1),ROI_temp2(temp4,2),ROI_temp2(temp4,3),10,'filled');

title(strcat('clust_',num2str(i),'_fmr1=',num2str(length(temp2)),"_avg",num2str(length(temp2)/11),'_hets=',num2str(length(temp5)),"_avg",num2str(length(temp5)/20),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/10)));


set(gcf, 'Position',  [100, 100, 450, 800])
view(-90,90);
end


%%% with the zbrain mask
%%% 11th loom
for i=[3 7 8 6]
    
    temp2=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,6},idx_temp2);
    temp5=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,6},idx_temp5);
    temp4=intersect(GoodClust_goodmembers_full(i).idx_thresh{1,6},idx_temp4);
    
figure; patch(brain3D,'EdgeColor','none','FaceAlpha',0.05);hold on;
% scatter3(ROI_temp2(GoodClust_goodmembers_full(i).idx,1),ROI_temp2(GoodClust_goodmembers_full(i).idx,2),ROI_temp2(GoodClust_goodmembers_full(i).idx,3),10,'filled');
% hold on;
scatter3(ROI_temp2(temp2,1),ROI_temp2(temp2,2),ROI_temp2(temp2,3),10,'filled');
hold on;
scatter3(ROI_temp2(temp5,1),ROI_temp2(temp5,2),ROI_temp2(temp5,3),10,'filled');
hold on;
scatter3(ROI_temp2(temp4,1),ROI_temp2(temp4,2),ROI_temp2(temp4,3),10,'filled');

title(strcat('clust_',num2str(i),'_fmr1=',num2str(length(temp2)),"_avg",num2str(length(temp2)/11),'_hets=',num2str(length(temp5)),"_avg",num2str(length(temp5)/20),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/10)));


set(gcf, 'Position',  [100, 100, 450, 800])
view(-90,90);
end