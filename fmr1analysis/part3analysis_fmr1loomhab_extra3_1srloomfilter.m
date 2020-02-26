%%% this script is based on the testing_deltaNcorr.m used for the
%%% habituation datasets.

%%% this script is to look at the ROIs that responded to the first loom.
%%% from there I can compare it to other looms and see the habituation but
%%% what I first will try is to look at the distribution of the max
%%% responses of the first 5 looms and the 11th. 

%%% I will do it in f20s 

load('s20_fmr1_loomhab_CN.mat','MatFiles','ZS_CN');


load('fmr1loomhab_BrainRegNclean.mat');


%%% to get the clusters from the Kmeans done to ALL the ROIs
load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN');

load('s20_good_idx_Fish.mat','idx_Fish');
idx_Fish_cat=categorical(idx_Fish);

load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');


load('s20_fmr1_loomhab_CN_part2_High_corr_Nb.mat','High_corr_Nb');

 %%%% i also need to generate the list of fish saved in the
 %%%% fmr1loomhab_lists.m file
%%% or load idx_temp2=fmr1; idx_temp4=controls and idx_temp5=hets
load('s20_fmr1_loomhab_CN_part3.mat','idx_temp2','idx_temp4','idx_temp5'); 


load('zbrain3D.mat');

%% generate a regressor for the first loom



Stimuli=zeros(1,100);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=60;
for i=1%:10
    Stimuli(i,(idxStart+(i-1)*60):(idxStart+(i-1)*60)+size(GCaMP6,1)-1)=GCaMP6/2;
end

figure;plot(Stimuli);
hold on;
plot(mean(ZS_CN(idx_rsq_cleaned,1:100)));

%% I will try with a linear regression to the first loom

    
ModelResults=[];
parfor i=1:size(ZS_CN,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(Stimuli',ZS_CN(i,1:100));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
     
end

rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them

figure; histogram(rsquare_loom);
std_rsq=std(rsquare_loom);
mean_rsq=mean(rsquare_loom);

idx_rsq1stloom2sd=find(rsquare_loom>mean_rsq+2*std_rsq & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
%%% 2SD seems to let a lot of crap still go through...


idx_rsq1stloom3sd=find(rsquare_loom>mean_rsq+3*std_rsq & rsquare_loom<1);
%%% 3SD is not bad at all...

figure; 
imagesc(ZS_CN(idx_rsq1stloom3sd,:), [-0.5 4]);colormap hot %%%to plot them in a raster plot

%%%% this filtering gets may of the interesting responses but when i plot
%%%% the location of the ROIs i see that they are still many outside the brain.  
figure;scatter(ROI_temp2(idx_rsq1stloom3sd,1),ROI_temp2(idx_rsq1stloom3sd,2),'filled');
hold on;
scatter(ROI_temp2(idx_rsq_cleaned,1),ROI_temp2(idx_rsq_cleaned,2),'filled');

figure;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.05);hold on;
set(gcf, 'Position',  [100, 100, 300, 600])
view(-90,90);
scatter3(ROI_temp2(idx_rsq1stloom3sd,1),ROI_temp2(idx_rsq1stloom3sd,2),ROI_temp2(idx_rsq1stloom3sd,3),10,'filled');

%%% what if I do 0.5? this on the other hand might be more restringent... i
%%% get far less ROIs than my idx_rsq_cleaned. lets see what happens
idx_rsq1stloom05=find(rsquare_loom>0.5 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; 
imagesc(ZS_CN(idx_rsq1stloom05,:), [-0.5 4]);colormap hot
figure;scatter(ROI_temp2(idx_rsq1stloom05,1),ROI_temp2(idx_rsq1stloom05,2),'filled');


%%%% its interesting because the fmr1 present a much stronger recovery
%%%% although the effect on the 2nd loom is not as clear. 
figure;
plot(mean(ZS_CN(intersect(idx_rsq1stloom3sd,idx_temp2),:)));
hold on;
plot(mean(ZS_CN(intersect(idx_rsq1stloom3sd,idx_temp5),:)));
hold on;
plot(mean(ZS_CN(intersect(idx_rsq1stloom3sd,idx_temp4),:)));
hold off;


%% to look at the distributions to the responses of the different looms but filtered after the first loom.



MaxPerFish_rsq=[];


fish=unique(idx_Fish);
list5=union(list1,list3);
edges=[0:0.25:15];
%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know why.
for f=1:length(unique(idx_Fish))
    tempfish=find(idx_Fish==fish(f));
      
      temp=intersect(idx_rsq1stloom3sd,tempfish);
      
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


%%% this part is to get  the means of the
%%% distribution to plot them. It seems that the fmr1 have stronger
%%% responses in the 2nd and the 11th loom!!!! although the hets look weird
%%% on the first loom... 

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

%%% and this to plot it in shades of grey and only first 3 looms and 11th

edges=[0.25:0.25:15];
group = fieldnames(MaxPerFish_rsq.MaxCountPerFish);
counter=1;

 figure;

   Whole_BinnedRespStrLooms_rsq_perFish=[];
   
for loom=[1 2 3 4 5 6];
subplot(1,6,counter);
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

%% locating the ROIs that respond more in 2nd and 11th loom

%%% I will do it with arbitrary thresholds

fish=unique(idx_Fish);
for f=1:length(unique(idx_Fish))
    tempfish=find(idx_Fish==fish(f));
      
      temp=intersect(idx_rsq1stloom3sd,tempfish);
      
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
    temp_std=std(idx_temp_res);
    temp_mean=mean(idx_temp_res);
    idx_thresh_resp=find(idx_temp_res>10);
    MaxPerFish_rsq.MaxIdxPerFish.(group).loom1.(strcat('idx_thresh_',num2str(fish(f))))=temp(idx_thresh_resp);
    %%% with mean+2SD i actually have more ROIs in controls... also for >10
    %%% but not as dramatic.
    
    idx_temp_res2=max(ZS_CN(temp,110:130),[],2);
    temp_std=std(idx_temp_res2);
    temp_mean=mean(idx_temp_res2);
    idx_thresh_resp=find(idx_temp_res2>3.25);
    MaxPerFish_rsq.MaxIdxPerFish.(group).loom2.(strcat('idx_thresh_',num2str(fish(f))))=temp(idx_thresh_resp);
    
    idx_temp_res2=max(ZS_CN(temp,145:165),[],2);
    temp_std=std(idx_temp_res2);
    temp_mean=mean(idx_temp_res2);
    idx_thresh_resp=find(idx_temp_res2>3.25);
    MaxPerFish_rsq.MaxIdxPerFish.(group).loom3.(strcat('idx_thresh_',num2str(fish(f))))=temp(idx_thresh_resp);
    
    idx_temp_res2=max(ZS_CN(temp,185:205),[],2);%Loom 4
    temp_std=std(idx_temp_res2);
    temp_mean=mean(idx_temp_res2);
    idx_thresh_resp=find(idx_temp_res2>3.25);
    MaxPerFish_rsq.MaxIdxPerFish.(group).loom4.(strcat('idx_thresh_',num2str(fish(f))))=temp(idx_thresh_resp);
    
    idx_temp_res2=max(ZS_CN(temp,220:240),[],2);%Loom 5
    temp_std=std(idx_temp_res2);
    temp_mean=mean(idx_temp_res2);
    idx_thresh_resp=find(idx_temp_res2>3.25);
    MaxPerFish_rsq.MaxIdxPerFish.(group).loom5.(strcat('idx_thresh_',num2str(fish(f))))=temp(idx_thresh_resp);
    
    idx_temp_res2=max(ZS_CN(temp,500:550),[],2);%Loom 11
    temp_std=std(idx_temp_res2);
    temp_mean=mean(idx_temp_res2);
    idx_thresh_resp=find(idx_temp_res2>5.5);
    MaxPerFish_rsq.MaxIdxPerFish.(group).loom6.(strcat('idx_thresh_',num2str(fish(f))))=temp(idx_thresh_resp);
    
end

for f=1:length(unique(idx_Fish))

 if ismember(fish(f),list2)
        group='fmr1';
        
    elseif ismember(fish(f),list4)
        group='control';
        
    elseif ismember(fish(f),list5)
        group='hets';
        
    else 
 end   
  
for loom=1:6    

MaxPerFish_rsq.MaxIdxPerFish.(group).(strcat('loom',num2str(loom))).idx_thresh_all=[];

end
end


for f=1:length(unique(idx_Fish))

 if ismember(fish(f),list2)
        group='fmr1';
        
    elseif ismember(fish(f),list4)
        group='control';
        
    elseif ismember(fish(f),list5)
        group='hets';
        
    else 
 end   
  
for loom=1:6    
temp_idx=MaxPerFish_rsq.MaxIdxPerFish.(group).(strcat('loom',num2str(loom))).(strcat('idx_thresh_',num2str(fish(f))));
MaxPerFish_rsq.MaxIdxPerFish.(group).(strcat('loom',num2str(loom))).idx_thresh_all=vertcat(MaxPerFish_rsq.MaxIdxPerFish.(group).(strcat('loom',num2str(loom))).idx_thresh_all,temp_idx);

end
end

%%%% IMPORTANT NOTE: I just realised that almost all ROIs of highly
%%%% respondent neurons in both the 2nd and the 11th looms belong to 1 fmr1
%%%% fish... is not the same one for 2nd than 11th. I think this is a big
%%%% problem... 

%%% to plot it in the brain



for loom=1:6
    
   temp2=MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
figure;
scatter(ROI_temp2(idx_rsq1stloom3sd,2),ROI_temp2(idx_rsq1stloom3sd,1),10,'filled');
hold on;
scatter(ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all,2),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all,1),10,'filled');
hold on;
scatter(ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all,2),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all,1),10,'filled');
hold on;
scatter(ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all,2),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all,1),10,'filled');

title(strcat('loom',num2str(loom),'_fmr1=',num2str(length(temp2)),"_avg",num2str(length(temp2)/11),'_hets=',num2str(length(temp5)),"_avg",num2str(length(temp5)/20),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/10)));

end

%%%%% with the zbrain mask

for loom=1:6
    
   temp2=MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
    
figure; patch(brain3D,'EdgeColor','none','FaceAlpha',0.05);hold on;

% scatter3(ROI_temp2(idx_rsq1stloom3sd,1),ROI_temp2(idx_rsq1stloom3sd,2),,ROI_temp2(idx_rsq1stloom3sd,3)10,'filled');
% hold on;
scatter3(ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all,1),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all,2),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all,3),10,'filled');
hold on;
% scatter3(ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all,1),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all,2),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all,3),10,'filled');
% hold on;
scatter3(ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all,1),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all,2),ROI_temp2(MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all,3),10,'filled');

title(strcat('loom',num2str(loom),'_fmr1=',num2str(length(temp2)),"_avg",num2str(length(temp2)/11),'_hets=',num2str(length(temp5)),"_avg",num2str(length(temp5)/20),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/10)));

set(gcf, 'Position',  [100, 100, 300, 600])
view(-90,90);
end



%%% to plot their means
for loom=1:6
    
   temp2=MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
figure;
plot(mean(ZS_CN(idx_rsq1stloom3sd,:)));
hold on;
plot(mean(ZS_CN(temp2,:)));
hold on;
plot(mean(ZS_CN(temp5,:)));
hold on;
plot(mean(ZS_CN(temp4,:)));

title(strcat('loom',num2str(loom),'_fmr1=',num2str(length(temp2)),"_avg",num2str(length(temp2)/11),'_hets=',num2str(length(temp5)),"_avg",num2str(length(temp5)/20),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/10)));

end






%% to find out if the highly respondent cells in the 2nd loom are the same ones than in the 11th... 

loom=6;
loom11_temp2=MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
loom11_temp5=MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
loom11_temp4=MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all;

StrongResp11_n_others=struct;
for loom=2:5
    
   temp2=MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
   StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all=intersect(loom11_temp2,temp2);
   StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all=intersect(loom11_temp5,temp5);
   StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all=intersect(loom11_temp4,temp4);

   
end

for loom=2:5
    
   temp2=StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
figure;
scatter(ROI_temp2(idx_rsq1stloom3sd,2),ROI_temp2(idx_rsq1stloom3sd,1),10,'filled');
hold on;
scatter(ROI_temp2(temp2,2),ROI_temp2(temp2,1),10,'filled');
hold on;
scatter(ROI_temp2(temp5,2),ROI_temp2(temp5,1),10,'filled');
hold on;
scatter(ROI_temp2(temp4,2),ROI_temp2(temp4,1),10,'filled');

title(strcat('loom',num2str(loom),'_fmr1=',num2str(length(temp2)),"_ratio",num2str(length(temp2)/length(loom11_temp2)),'_ratio=',num2str(length(temp5)),"_ratio",num2str(length(temp5)/length(loom11_temp5)),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/length(loom11_temp4))));

end

%%% with zbrain mask
for loom=2:5
    
   temp2=StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all;
    
figure; patch(brain3D,'EdgeColor','none','FaceAlpha',0.05);hold on;


scatter3(ROI_temp2(temp2,1),ROI_temp2(temp2,2),ROI_temp2(temp2,3),10,'filled');
hold on;
scatter3(ROI_temp2(temp5,1),ROI_temp2(temp5,2),ROI_temp2(temp5,3),10,'filled');
hold on;
scatter3(ROI_temp2(temp4,1),ROI_temp2(temp4,2),ROI_temp2(temp4,3),10,'filled');

title(strcat('loom',num2str(loom),'_fmr1=',num2str(length(temp2)),"_ratio",num2str(length(temp2)/length(loom11_temp2)),'_ratio=',num2str(length(temp5)),"_ratio",num2str(length(temp5)/length(loom11_temp5)),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/length(loom11_temp4))));

set(gcf, 'Position',  [100, 100, 300, 600])
view(-90,90);
end


%%% to plot their means
for loom=2:5
    
   temp2=StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
   
figure;
plot(mean(ZS_CN(idx_rsq1stloom3sd,:)));
hold on;
plot(mean(ZS_CN(temp2,:)));
hold on;
plot(mean(ZS_CN(temp5,:)));
hold on;
plot(mean(ZS_CN(temp4,:)));

title(strcat('loom',num2str(loom),'_fmr1=',num2str(length(temp2)),"_avg",num2str(length(temp2)/11),'_hets=',num2str(length(temp5)),"_avg",num2str(length(temp5)/20),'_wt=',num2str(length(temp4)),"_avg",num2str(length(temp4)/10)));

end


%% to see which cluster they belong to


%%% for all the strong respondent ROIs

for loom=1:6
    
   temp2=MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID=idxKmeans_ZS_CN(temp2);
MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID=idxKmeans_ZS_CN(temp5);
MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID=idxKmeans_ZS_CN(temp4);

end

for loom=1:6

    temp2=unique(MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID);
    temp5=unique(MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID);
    temp4=unique(MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID);

    
        for i=1:length(temp2)
        MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count(temp2(i))=length(find(temp2(i)==MaxPerFish_rsq.MaxIdxPerFish.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID));
        end

        for i=1:length(temp5)
        MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count(temp5(i))=length(find(temp5(i)==MaxPerFish_rsq.MaxIdxPerFish.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID));
        end
        
        for i=1:length(temp4)
        MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count(temp4(i))=length(find(temp4(i)==MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID));
        end
        
end


for loom=1:6
figure;bar(MaxPerFish_rsq.MaxIdxPerFish.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count);
end


%%%% so it seems that must ROIs belong to clusters 36(slopehab like)
%%%% and 40 (sharp fasthab like). so thats interesting. 
%%%% but for fmr1 the 11th loom also has from cluster 13 (the weird one).



%%% to check which fish has more ROIs in the odd cluster (13) and see if
%%% its fmr1 fish...
figure;histogram(idx_Fish_cat(find(idxKmeans_ZS_CN==13)));

%%%% it seems that most of it belongs to fish 201811081 (more than 2500
%%%% ROIs). the other fish with most is 201811087 (almost 1500). the rest
%%%% of the fish have around 200 ROIs.

%%%% 201811087 is an fmr1 fish and is the one that brings the most ROIs to
%%%% the loom responses to the 11th loom for that genotype.

%%% to see the responses of the neurons of that fish from cluster 13. 
temp_idx1=find(idxKmeans_ZS_CN==13);
temp_idx2=intersect(temp_idx1,find(idx_Fish==201811087));
figure;plot(mean(ZS_CN(temp_idx2,:)));


%%%% 201811081 is a control fish!!! and only brings 11 ROIs to the 11th
%%%% respondent group. probably because it didnt responde much at the
%%%% begining... 
temp_idx1=find(idxKmeans_ZS_CN==13);
temp_idx2=intersect(temp_idx1,find(idx_Fish==201811081));
figure;plot(mean(ZS_CN(temp_idx2,:)));

%%% IMPORTANT NOTE: maybe what I should do is to get the top ~2%(2SD) ROIs of
%%% each fish for each loom and compare them from that... Or a specific
%%% number of ROIs per fish... not sure how to proceed. 



%%% for all the overlaping strong respondent ROIs

for loom=2:5
    
   temp2=StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp5=StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all;
   temp4=StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all;
   
StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID=idxKmeans_ZS_CN(temp2);
StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID=idxKmeans_ZS_CN(temp5);
StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID=idxKmeans_ZS_CN(temp4);

end

for loom=2:5

    temp2=unique(StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID);
    temp5=unique(StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID);
    temp4=unique(StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID);

    
        for i=1:length(temp2)
        StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count(temp2(i))=length(find(temp2(i)==StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID));
        end

        for i=1:length(temp5)
       StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count(temp5(i))=length(find(temp5(i)==StrongResp11_n_others.hets.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID));
        end
        
        for i=1:length(temp4)
        StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count(temp4(i))=length(find(temp4(i)==StrongResp11_n_others.control.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID));
        end
        
end


for loom=2:5
figure;bar(StrongResp11_n_others.fmr1.(strcat('loom',num2str(loom))).idx_thresh_all_KmeansID_count);
end

%%%% the results for the overlap betwen 2nd loom and 11th loom are
%%%% interesting. most of the Overlap is in the OT although a bit is in Dm.
%%%% and also interestinglgy most neurons belong to cluster 36(slopehab)
%%%% and then to cluster 40 (sharp fasthab) for fmr1 fish and this is inverted for controls. 



