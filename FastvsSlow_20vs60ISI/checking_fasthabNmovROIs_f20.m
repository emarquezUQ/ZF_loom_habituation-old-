

%%% this script is to check again on the tail movements. I will test it
%%% with f20 and/or s20 because the 60ISI movies had some issues picking up
%%% the ROIs with movements. 


load('All_More_BrainReg2.mat')
load('Zbrain_brainMask2D.mat')

f20_idx=load('f20_cleaned_idxs.mat');
%s20_idx=load('s20_cleaned_idxs.mat');


load('ZS_N_idx_Fish_all.mat','ZS_f20','idx_Fish_f20');


%%%% to check the movement responses per fish. 
counter=1;
figure;
for i=unique(idx_Fish_f20)'
    subplot(3,4,counter);
    temp_idx=find(idx_Fish_f20==i);
    temp_idx=intersect(temp_idx,f20_idx.idx_rsq_Mov_cleaned);
    imagesc(ZS_f20(temp_idx,:));
counter=counter+1;
    
end

%% first checking which green neurons overlap with the movement neurons. 

fasthab_idx=f20_idx.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned;

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(ROI_temp2.f20(fasthab_idx,1),ROI_temp2.f20(fasthab_idx,2),10,'filled')
view(-90,90)



fastNmov_idx=intersect(fasthab_idx,f20_idx.idx_rsq_Mov_cleaned);

%%% how many and wich fish are included? the usual ones: 1, 13, 34 and 44.
unique(idx_Fish_f20(fastNmov_idx))
unique(idx_Fish_f20)

%%% where are these neurons? a bit all over. 
figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(ROI_temp2.f20(fastNmov_idx,1),ROI_temp2.f20(fastNmov_idx,2),10,'filled')
view(-90,90)


%% a correlation
%%%% now, a correlation of each fish fasthab ROIs to its movements to see if i can find more neurons related

%%%% I will not count the first loom cause it could influence too much on
%%%% the results

%%% NOTE: my r2 filter of 0.2 to filter the mov ROIs is more or less equivalent to a corr coef of 0.447 


%%%%%% this bit was taking too long. I might tryagain later... 

% movCorr=zeros(size(idx_Fish_f20));
% for i=1:length(ZS_f20)
%     
%     temp_f=idx_Fish_f20(i);
%     temp_idx=find(idx_Fish_f20==temp_f);
%     temp_idx=intersect(temp_idx,f20_idx.idx_rsq_Mov_cleaned);
%     temp_mov=mean(ZS_f20(temp_idx,:));
%     
%     temp_corr=corrcoef(ZS_f20(i,120:end),temp_mov(1,120:end));
%     temp_corr=temp_corr(1,2);
%     
%    movCorr(i)=temp_corr; 
%     
% end
% 
% 
% idx_movCorr=find(movCorr>0.5); %%% i first tested it with 0.5
% 

% figure;
% plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
% hold on;
% scatter(ROI_temp2.f20(idx_movCorr,1),ROI_temp2.f20(idx_movCorr,2),10,'filled')
% view(-90,90)
% 
% 
% fastNmovCorr_idx=intersect(fasthab_idx,idx_movCorr);


%%

%%% first checking normality, to see if I use Pearson of Spearman. which i
%%% think are both wrong anyway because they are not independent... 
normTest=zeros(length(idx_Fish_f20(fasthab_idx)),3);
for i=1:length(idx_Fish_f20(fasthab_idx))
    
    temp_f=idx_Fish_f20(fasthab_idx(i));
    temp_idx=find(idx_Fish_f20==temp_f);
        
    tempKS=kstest(ZS_f20(fasthab_idx(i),120:end));
    %tempAD=adtest(ZS_f20(fasthab_idx(i),120:end));
    tempJB=jbtest(ZS_f20(fasthab_idx(i),120:end));
        
    normTest(i,1)=tempKS; 
    %normTest(i,2)=tempAD; 
    normTest(i,2)=tempJB;
       
end

%%% so it seems that most of them are not normal so i will use spearman 

fastNmovCorr=zeros(length(idx_Fish_f20(fasthab_idx)),2);
for i=1:length(idx_Fish_f20(fasthab_idx))
    
    temp_f=idx_Fish_f20(fasthab_idx(i));
    temp_idx=find(idx_Fish_f20==temp_f);
    temp_idx=intersect(temp_idx,f20_idx.idx_rsq_Mov_cleaned);
    temp_mov=mean(ZS_f20(temp_idx,:));
    
    [temp_corr temp_p]=corr(ZS_f20(fasthab_idx(i),120:end)',temp_mov(1,120:end)','Type','Spearman');
    temp_corr=temp_corr(1,1);
    temp_p=temp_p(1,1);
    
   fastNmovCorr(i,1)=temp_corr; 
    fastNmovCorr(i,2)=temp_p;  
end

figure;histogram(fastNmovCorr(:,1));

h = kstest(fastNmovCorr(:,1))
h = adtest(fastNmovCorr(:,1))
h = jbtest(fastNmovCorr(:,1))

length(find(isnan(fastNmovCorr(:,1))));%%% there are some NANs cause i had fish without movements. 
nanmean(fastNmovCorr(:,1));
nanstd(fastNmovCorr(:,1));

threshold=nanmean(fastNmovCorr(:,1))+nanstd(fastNmovCorr(:,1))*1; %%% with Pearson 0.5561 with 2SD and 0.3650 with 1SD; with spearman 0.4610 for 2SD and 0.3066 fro 1SD

%%% to selected based on corr coef value
idx_fastNmovCorr=find(fastNmovCorr(:,1)>threshold); %%% i first tested it with 0.5
idx_fastNmovCorr=fasthab_idx(idx_fastNmovCorr);

%%% to selected based on significance. it works pretty bad because most of
%%% them are significant
% idx_fastNmovCorr=find(fastNmovCorr(:,2)<0.05); 
% idx_fastNmovCorr=fasthab_idx(idx_fastNmovCorr);


%%% where are these neurons? is looking good...
figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(ROI_temp2.f20(idx_fastNmovCorr,1),ROI_temp2.f20(idx_fastNmovCorr,2),10,'filled')
view(-90,90)

% figure;
% scatter3(ROI_temp2.f20(idx_fastNmovCorr,1),ROI_temp2.f20(idx_fastNmovCorr,2),ROI_temp2.f20(idx_fastNmovCorr,3),10,'filled')
% view(-90,90)


%%%% comparing them with the total fast hab
figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(ROI_temp2.f20(fasthab_idx,1),ROI_temp2.f20(fasthab_idx,2),10,'filled')
hold on;
scatter(ROI_temp2.f20(idx_fastNmovCorr,1),ROI_temp2.f20(idx_fastNmovCorr,2),10,'filled')
view(-90,90)



%%% how many and wich fish are included? 7-10 of them... not
%%% bad!

unique(idx_Fish_f20(idx_fastNmovCorr))
%unique(idx_Fish_f20)

%%%% how many are the same that are already overlaping? not all... 
length(intersect(idx_fastNmovCorr,fastNmov_idx))



%%%% Ok... but how do the actual resposnes look like?

%%%% to check the movement responses per fish. 
counter=1;
figure;
for i=unique(idx_Fish_f20)'
    subplot(3,4,counter);
    temp_idx=find(idx_Fish_f20==i);
    temp_idx=intersect(temp_idx,idx_fastNmovCorr);
    imagesc(ZS_f20(temp_idx,:));
counter=counter+1;
    
end

%%%% to check the movement responses per fish in the brain. 
counter=1;
figure;
for i=unique(idx_Fish_f20)'
    subplot(3,4,counter);
    temp_idx=find(idx_Fish_f20==i);
    temp_idx=intersect(temp_idx,idx_fastNmovCorr);
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROI_temp2.f20(temp_idx,1),ROI_temp2.f20(temp_idx,2),10,'filled')
    view(-90,90)
       
counter=counter+1;
    
end




%%% the means together with the means of the Mov ROIs. and the nonhab to
%%% see when do loom happens

nonhabMean=mean(ZS_f20(f20_idx.clust_f20_CL4_cleaned.clust_f20_CL4_1_cleaned,:));

counter=1;
figure;
for i=unique(idx_Fish_f20)'
    subplot(3,4,counter);
    f_idx=find(idx_Fish_f20==i);
    temp_idx=intersect(f_idx,idx_fastNmovCorr);
    plot(mean(ZS_f20(temp_idx,:)));
    hold on;
    plot(mean(ZS_f20(intersect(f_idx,f20_idx.idx_rsq_Mov_cleaned),:)));
    hold on;
    plot(nonhabMean);
counter=counter+1;
    
end



%%%% to see how many ROIs are in each brain region. 

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

fasthab_mov_brain=[];
 for brain=1:length(RegionList)   
     
      idx_brain_temp=PerBrainRegions.f20.(RegionList{brain}).idx;
      temp_idx=intersect(idx_fastNmovCorr,idx_brain_temp);
      
      if ~isempty(temp_idx)
      fasthab_mov_brain(brain)=size(temp_idx,1);
      
      else       
        fasthab_mov_brain(brain)=0;  
      end
 end
 
 figure;bar(fasthab_mov_brain);
 
 
 %%% and the proportions compared to the fasthab with in each brain region.
 fasthab_mov_brain_prop=[];
 for brain=1:length(RegionList)   
     
      idx_brain_temp=PerBrainRegions.f20.(RegionList{brain}).idx;
      temp_idx1=intersect(idx_fastNmovCorr,idx_brain_temp);
      temp_idx2=intersect(fasthab_idx,idx_brain_temp);
      
      if ~isempty(temp_idx1)
      fasthab_mov_brain_prop(brain)=size(temp_idx1,1)/size(temp_idx2,1);
      
      else       
        fasthab_mov_brain_prop(brain)=0;  
      end
 end
 
 figure;bar(fasthab_mov_brain_prop);
 
 %%% proportion compared to the whole brain fasthab ROIs... not very useful
%  fasthab_mov_brain_prop=[];
%  for brain=1:length(RegionList)   
%      
%       idx_brain_temp=PerBrainRegions.f20.(RegionList{brain}).idx;
%       temp_idx1=intersect(idx_fastNmovCorr,idx_brain_temp);     
%       
%       if ~isempty(temp_idx1)
%       fasthab_mov_brain_prop(brain)=size(temp_idx1,1)/size(fasthab_idx,1);
%       
%       else       
%         fasthab_mov_brain_prop(brain)=0;  
%       end
%  end
%  
%  figure;bar(fasthab_mov_brain_prop);
 
 
 
 
 %%%% to plot them in the brain colored by brain region. 
 counter=1;
figure;
for brain=1:length(RegionList)  
     idx_brain_temp=PerBrainRegions.f20.(RegionList{brain}).idx;
      temp_idx=intersect(idx_fastNmovCorr,idx_brain_temp);
    scatter(ROI_temp2.f20(temp_idx,1),ROI_temp2.f20(temp_idx,2),10,'filled')
    view(-90,90)
       hold on;
counter=counter+1;
    
end
 


%%%% getting the ROIs locations for Unity
CSV_temp=ROI_temp2.f20(idx_fastNmovCorr,:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_fastNmovCorr_f20_','.csv');

    csvwrite(filename,CSV_temp);



%%%% I will try this later. 

%%%%% ok, now it seems that a subset of fasthab neurons do correlate with
%%%%% behavioural outputs. what happens in the matrices when those movents
%%%%% coincide with a loom? lest try to have a look
load('graph_analysis_loomhab3.mat')
load('graph_loomNdataset2.mat')

%%%% i should try per individual fish... but many dont have all the nodes.
%%%% i will try first if i can detect when the movements happen

datasets=fieldnames(Data_corrMat2);

 data=1
    
 datatemp=datasets{data};

 fish=fieldnames(Data_corrMat2.(datatemp));
 
  mean_corr_mat=[];
 for f=1:length(fish)-1

 for i=1:31
 temp_mat=Data_corrMat2.(datatemp).(fish{f}).loomsR{1,i}(keep,keep);
 temp_mat(find(isnan(temp_mat)))=0;
 mean_corr_mat(f,i)=mean(mean(temp_mat));
 
 end
 end
 
 figure;
 for i=1:11
 plot(mean_corr_mat(i,:));
 hold on;
 end
 
 figure;
 plot(mean(mean_corr_mat));

 