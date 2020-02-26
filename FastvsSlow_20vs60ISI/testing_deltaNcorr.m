

%%% this script is to look at the difference in response to the first loom
%%% and to the 10th loom to make a heat map of this change. I will select 
%%% ROIs that responded to the first loom and then look how they changed.
%%% I will select this ROIs with dF/F to a 2SD response, a simpble
%%% correlation and a linear regression. They will give similar results but
%%% is to see which one looks better. 

%%% I will do it in f20s 

ZS_f20=load('f20_CN_r2050_CL5_extra.mat','ZS_CN');
 

load('All_More_BrainReg2.mat');


f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');


%%

%%% this is if I want to get first the ROIs that reacted when the first
%%% loom presentation happen based on the strenght of the response
%%% filtering by 2 units of the z-score. then I get the delta of the first
%%% response and substract it to the delta of the 10th response
idx_resp=[];
for i=1:size(ZS_f20.ZS_CN,1)

%if  max(ZS_f20.ZS_CN(i,60:80))>2   
idx_resp(i)=max(ZS_f20.ZS_CN(i,60:80))>2;
%elseif min(ZS_f20.ZS_CN(i,60:80))<-2 
%idx_resp(i)=min(ZS_f20.ZS_CN(i,60:80))<-2;
%end

delta_1(i)=max(ZS_f20.ZS_CN(i,60:80))-mean(ZS_f20.ZS_CN(i,10:59));
delta_2(i)=max(ZS_f20.ZS_CN(i,105:140))-mean(ZS_f20.ZS_CN(i,10:59));
delta_10(i)=max(ZS_f20.ZS_CN(i,420:440))-mean(ZS_f20.ZS_CN(i,10:59));

deltahab2(i)=delta_1(i)-delta_2(i);
deltahab10(i)=delta_1(i)-delta_10(i);
 
ratiohab2(i)=delta_2(i)/delta_1(i);
ratiohab10(i)=delta_10(i)/delta_1(i);

end

idx_resp=find(idx_resp);

figure;imagesc(ZS_f20.ZS_CN(idx_resp,:)); colormap('hot');
figure;histogram(deltahab2);
figure;histogram(deltahab10);

%% 
%%% then I get the SD of the difference between 1st and 10th to select the 
%%% ROIs that had above or bellow a 2 SD difference

std_delta=std(deltahab10);

idx_delta_pos=find(deltahab10>std_delta*2);
idx_delta_neg=find(deltahab10<-std_delta*2);
idx_delta=union(idx_delta_neg,idx_delta_pos);

histogram(deltahab10(idx_delta));

%%% this is to look at the negative ones. they look like noise
figure;imagesc(ZS_f20.ZS_CN(idx_delta_neg,:));  %%%% to me they look more like noise... 
figure;plot(mean(ZS_f20.ZS_CN(idx_delta_neg,:)));  %%%% to me they look more like noise... 


figure;imagesc(ZS_f20.ZS_CN(idx_delta_pos,:));  
figure;plot(mean(ZS_f20.ZS_CN(idx_delta_pos,:)));  


%%% so i will only use the positive ones

deltahab_2std=deltahab10(idx_delta_pos);
Qs = quantile(deltahab_2std,[0.025 0.25 0.50 0.75 0.975]);

figure;
scatter(ROI_temp2.f20(idx_delta_pos,1),ROI_temp2.f20(idx_delta_pos,2),10,deltahab_2std,'filled');colormap('jet');colorbar; caxis([Qs(1) Qs(5)]);
figure;
scatter3(ROI_temp2.f20(idx_delta_pos,1),ROI_temp2.f20(idx_delta_pos,2),ROI_temp2.f20(idx_delta_pos,3),10,deltahab_2std,'filled');colormap('jet');colorbar; caxis([Qs(1) Qs(5)]);
 
%%%%% it looks ok... but a bit noisy. 

%%

%%%% now i will filter the ROIs with a correlation

%%

%%% now i will try to do the same but with a correlation

Stimuli=zeros(1,100);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=60;
for i=1%:10
    Stimuli(i,(idxStart+(i-1)*60):(idxStart+(i-1)*60)+size(GCaMP6,1)-1)=GCaMP6;
end

figure;plot(Stimuli);

for i=1:size(ZS_f20.ZS_CN,1)

temp_corr=corrcoef(ZS_f20.ZS_CN(i,1:100),Stimuli);
Corr_f20(i)=temp_corr(1,2);

end

figure;histogram(Corr_f20);

%idx_corr05=find(Corr_f20>0.5);


std_corr=std(Corr_f20);

idx_corr_pos=find(Corr_f20>std_corr*2);
idx_corr_neg=find(Corr_f20<-std_corr*2);
idx_corr=union(idx_corr_neg,idx_corr_pos);

histogram(Corr_f20(idx_corr));
figure;imagesc(ZS_f20.ZS_CN(idx_corr,:)); 
figure;plot(mean(ZS_f20.ZS_CN(idx_corr,:))); 

figure;imagesc(ZS_f20.ZS_CN(idx_corr_pos,:)); 
figure;plot(mean(ZS_f20.ZS_CN(idx_corr_pos,:))); 

figure;imagesc(ZS_f20.ZS_CN(idx_corr_neg,:));   %%% this ones look more real... but their locations not that much
figure;plot(mean(ZS_f20.ZS_CN(idx_corr_neg,:))); %%% this ones look more real... but their locations not that much

figure;
scatter3(ROI_temp2.f20(idx_corr_neg,1),ROI_temp2.f20(idx_corr_neg,2),ROI_temp2.f20(idx_corr_neg,3),10,idx_corr_neg,'filled');colormap('jet');colorbar; caxis([Qs_Corr(1) Qs_Corr(5)]);


Corr_2std=Corr_f20(idx_corr);
Qs_Corr = quantile(Corr_2std,[0.025 0.25 0.50 0.75 0.975]);

figure;
scatter(ROI_temp2.f20(idx_corr,1),ROI_temp2.f20(idx_corr,2),10,Corr_2std,'filled');colormap('jet');colorbar; caxis([Qs_Corr(1) Qs_Corr(5)]);
figure;
scatter3(ROI_temp2.f20(idx_corr,1),ROI_temp2.f20(idx_corr,2),ROI_temp2.f20(idx_corr,3),10,Corr_2std,'filled');colormap('jet');colorbar; caxis([Qs_Corr(1) Qs_Corr(5)]);

%%

%%% if i use both criterions together? it looks cleaner. 

idx_DelCorr=intersect(idx_delta,idx_corr);

figure;histogram(Corr_f20(idx_DelCorr));
figure;imagesc(ZS_f20.ZS_CN(idx_DelCorr,:));   
figure;plot(mean(ZS_f20.ZS_CN(idx_DelCorr,:)));  



Corr_2std2=Corr_f20(idx_DelCorr);
Qs_Corr2 = quantile(Corr_2std2,[0.025 0.25 0.50 0.75 0.975]);

figure;
scatter(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),10,Corr_2std2,'filled');colormap('jet');colorbar; caxis([Qs_Corr2(1) Qs_Corr2(5)]);
figure;
scatter3(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),ROI_temp2.f20(idx_DelCorr,3),10,Corr_2std2,'filled');colormap('jet');colorbar; caxis([Qs_Corr2(1) Qs_Corr2(5)]);

%%% ploting the heatmap with the deltas. it doesnt look bad. 
delNcorr=deltahab10(idx_DelCorr);


figure;
scatter(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),10,delNcorr,'filled');colormap('hot');colorbar; caxis([Qs(1) Qs(5)]);
figure;
scatter3(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),ROI_temp2.f20(idx_DelCorr,3),10,delNcorr,'filled');colormap('hot');colorbar; caxis([Qs(1) Qs(5)]);

%%% now trying to do a figure for 1st and for 10th loom. 
del_1_Ncorr=delta_1(idx_DelCorr);
Qs_d1 = quantile(del_1_Ncorr,[0.025 0.25 0.50 0.75 0.975]);

del_2_Ncorr=delta_10(idx_DelCorr);
Qs_d10 = quantile(del_2_Ncorr,[0.025 0.25 0.50 0.75 0.975]);



figure;
subplot(1,3,1)
scatter(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),10,del_1_Ncorr,'filled');colormap('hot');colorbar; caxis([Qs_d1(1) Qs_d1(5)]);
subplot(1,3,2)
scatter(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),10,del_2_Ncorr,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d10(5)]);
subplot(1,3,3)
scatter(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),10,delNcorr,'filled');colormap('hot');colorbar; caxis([Qs(1) Qs(5)]);


%%

%%%% now I will look at the correlation to the same regressor at the 10th
%%%% loom... I am not sure this is the best aproach cause a different
%%%% response would have low correlation but could be higher in
%%%% intensity...

for i=1:size(ZS_f20.ZS_CN,1)

temp_corr=corrcoef(ZS_f20.ZS_CN(i,420:440),Stimuli(60:80));
Corr_f20_10th(i)=temp_corr(1,2);

end

figure;histogram(Corr_f20_10th);


Corr_2std2_10th=Corr_f20_10th(idx_DelCorr);
Qs_Corr3 = quantile(Corr_2std2_10th,[0.025 0.25 0.50 0.75 0.975]);

figure;histogram(Corr_2std2_10th);

figure;
scatter(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),10,Corr_2std2_10th,'filled');colormap('jet');colorbar; caxis([Qs_Corr3(1) Qs_Corr3(5)]);
figure;
scatter3(ROI_temp2.f20(idx_DelCorr,1),ROI_temp2.f20(idx_DelCorr,2),ROI_temp2.f20(idx_DelCorr,3),10,Corr_2std2_10th,'filled');colormap('jet');colorbar; caxis([Qs_Corr3(1) Qs_Corr3(5)]);

%%%%%  not very informative. i think... 

%%
%%% finally, I will try with a linear regression to the first loom


    
ModelResults=[];
parfor i=1:size(ZS_f20.ZS_CN,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(Stimuli',ZS_f20.ZS_CN(i,1:100));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
     
end

rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them

figure; histogram(rsquare_loom);
std_rsq=std(rsquare_loom);

idx_rsq=find(rsquare_loom>2*std_rsq & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; 
imagesc(ZS_f20.ZS_CN(idx_rsq,:), [-0.5 4]);colormap hot %%%to plot them in a raster plot

%%%% this filtering gets may of the interesting responses but when i plot
%%%% the location of the ROIs i see that they are still many outside the brain.  

%%% what if I do 0.5?
idx_rsq=find(rsquare_loom>0.5 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; 
imagesc(ZS_f20.ZS_CN(idx_rsq,:), [-0.5 4]);colormap hot


rsq2d1_2=deltahab2(idx_rsq);
Qs_rsq_2 = quantile(rsq2d1_2,[0.025 0.25 0.50 0.75 0.975]);
figure;histogram(rsq2d1_2);
temp_std=std(rsq2d1_2)
figure;histogram(rsq2d1_2(find(rsq2d1_2<-2*temp_std)))
figure;plot(mean(ZS_f20.ZS_CN(idx_rsq(find(rsq2d1_2<-2*temp_std)),:))); %%% in this case it looks like noise
figure;imagesc(ZS_f20.ZS_CN(idx_rsq(find(rsq2d1_2<-2*temp_std)),:)); 


rsq2d1_10=deltahab10(idx_rsq);
Qs_rsq_10 = quantile(rsq2d1_10,[0.025 0.25 0.50 0.75 0.975]);
temp_std=std(rsq2d1_10)
figure;histogram(rsq2d1_10(find(rsq2d1_10<-2*temp_std)))
figure;plot(mean(ZS_f20.ZS_CN(idx_rsq(find(rsq2d1_10<-2*temp_std)),:))); %%% this is weird. I find some negative deltas (more response in the tenth loom) but when I plot them it doenst seem so 


rsq2r1_2=ratiohab2(idx_rsq);
Qs_rsq_ratio_2 = quantile(rsq2r1_2,[0.025 0.25 0.50 0.75 0.975]);
figure;histogram(rsq2r1_2); %%% i am geting some few but very weird results very far from 0-1. like -280.38 or 808.0434. they are most probably artifacts
max(rsq2r1_2)
high_idx=find(rsq2r1_2>1);
high=rsq2r1_2(high_idx);
figure;scatter(ROI_temp2.f20(idx_rsq(high_idx),1),ROI_temp2.f20(idx_rsq(high_idx),2),10,high,'filled');colormap('jet');colorbar; %caxis([Qs_rsq2(1) Qs_rsq2(5)]);



rsq2r1_10=ratiohab10(idx_rsq);
Qs_rsq_ratio_10 = quantile(rsq2r1_10,[0.025 0.25 0.50 0.75 0.975]);


del_1_Nrsq=delta_1(idx_rsq);
Qs_d1 = quantile(del_1_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_2_Nrsq=delta_2(idx_rsq);
Qs_d2 = quantile(del_2_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_10_Nrsq=delta_10(idx_rsq);
Qs_d10 = quantile(del_10_Nrsq,[0.025 0.25 0.50 0.75 0.975]);



figure;
subplot(1,5,1)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,del_1_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,2)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,del_2_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,3)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,del_10_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,4)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,rsq2r1_2,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_2(1) Qs_rsq_ratio_10(5)]);
subplot(1,5,5)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,rsq2r1_10,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_2(1) Qs_rsq_ratio_10(5)]);

%%
%%% delta of 2SD and rsq together. it looks clean. but this would be to
%%% look specifically for cells that had a big change 

idx_DelRsq=intersect(idx_delta,idx_rsq);


rsq2std2=deltahab10(idx_DelRsq);
Qs_rsq2 = quantile(rsq2std2,[0.025 0.25 0.50 0.75 0.975]);

del_1_Nrsq2=delta_1(idx_DelRsq);
Qs_d1_2 = quantile(del_1_Nrsq2,[0.025 0.25 0.50 0.75 0.975]);

del_2_Nrsq2=delta_10(idx_DelRsq);
Qs_d2_2 = quantile(del_2_Nrsq2,[0.025 0.25 0.50 0.75 0.975]);



figure;
subplot(1,3,1)
scatter(ROI_temp2.f20(idx_DelRsq,1),ROI_temp2.f20(idx_DelRsq,2),10,del_1_Nrsq2,'filled');colormap('hot');colorbar; caxis([Qs_d2_2(1) Qs_d1_2(5)]);
subplot(1,3,2)
scatter(ROI_temp2.f20(idx_DelRsq,1),ROI_temp2.f20(idx_DelRsq,2),10,del_2_Nrsq2,'filled');colormap('hot');colorbar; caxis([Qs_d2_2(1) Qs_d1_2(5)]);
subplot(1,3,3)
scatter(ROI_temp2.f20(idx_DelRsq,1),ROI_temp2.f20(idx_DelRsq,2),10,rsq2std2,'filled');colormap('jet');colorbar; caxis([Qs_rsq2(1) Qs_rsq2(5)]);

figure;
scatter3(ROI_temp2.f20(idx_DelRsq,1),ROI_temp2.f20(idx_DelRsq,2),ROI_temp2.f20(idx_DelRsq,3),10,rsq2std2,'filled');colormap('jet');colorbar; caxis([Qs_rsq2(1) Qs_rsq2(5)]);

%% to clean the ROIs outside the brain. of the linear regression filtering

%%% i first make a mask with all the brain regions
load('Zbrain_Masks.mat');


Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
 
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

figure;scatter(Zbrain_AllMask(:,1),Zbrain_AllMask(:,2),'.');


idx_brain=ismember(ROI_temp2.f20,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

%idxKmeans_final(~idx_brain)=0;  %%% this is how gilles had it but his only
                                 %%% work for indexes with all the ROIs. i need to adjust it.
 
idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

idx_rsq_1stLoom_cleaned=intersect(idx_brain2,idx_rsq);

rsq2d1_2=deltahab2(idx_rsq_1stLoom_cleaned);
Qs_rsq_2 = quantile(rsq2d1_2,[0.025 0.25 0.50 0.75 0.975]);

rsq2d1_10=deltahab10(idx_rsq_1stLoom_cleaned);
Qs_rsq_10 = quantile(rsq2d1_10,[0.025 0.25 0.50 0.75 0.975]);

rsq2r1_2=ratiohab2(idx_rsq_1stLoom_cleaned);
Qs_rsq_ratio_2 = quantile(rsq2r1_2,[0.025 0.25 0.50 0.75 0.975]);

rsq2r1_10=ratiohab10(idx_rsq_1stLoom_cleaned);
Qs_rsq_ratio_10 = quantile(rsq2r1_10,[0.025 0.25 0.50 0.75 0.975]);

del_1_Nrsq=delta_1(idx_rsq_1stLoom_cleaned);
Qs_d1 = quantile(del_1_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_2_Nrsq=delta_2(idx_rsq_1stLoom_cleaned);
Qs_d2 = quantile(del_2_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_10_Nrsq=delta_10(idx_rsq_1stLoom_cleaned);
Qs_d10 = quantile(del_10_Nrsq,[0.025 0.25 0.50 0.75 0.975]);



figure;
subplot(1,5,1)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,del_1_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,2)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,del_2_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,3)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,del_10_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,4)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_2,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_10(1) Qs_rsq_ratio_10(5)]);
subplot(1,5,5)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_10,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_10(1) Qs_rsq_ratio_10(5)]);


%%% for the colorbar of the substractions. I made a a colormap of 18 values
%%% (0:17) and upload it to unity. so I need that range for my colorbar. 

%c = jet(18);
%c = plasma(18);
%c = inferno(18);
c = viridis(18);


% filename=strcat('jet_colormap2.csv');
% filename=strcat('plasma_colormap.csv');
% filename=strcat('inferno_colormap.csv');
 filename=strcat('viridis_colormap.csv');
 
 csvwrite(filename,c);

figure;
 subplot(1,2,1)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_2,'filled');colormap(plasma);colorbar; caxis([0 1]);
subplot(1,2,2)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_10,'filled');colormap(plasma);colorbar; caxis([0 1]);


%% generate coordenates for Unity

%%%% I will put the coords of x,y and z. and then the rsq value and then
%%%% the deltas. 


%%% for the linear regression results

rsquare_loom_short=rsquare_loom(idx_rsq_1stLoom_cleaned);

    idx_temp2=idx_rsq_1stLoom_cleaned;
    
    CSV_temp=ROI_temp2.f20(idx_temp2,:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,5)=rsquare_loom_short;
    
%     CSV_temp(:,4)=del_1_Nrsq;
%     CSV_temp(:,4)=del_2_Nrsq;
%     CSV_temp(:,4)=del_10_Nrsq;
%     CSV_temp(:,4)=rsq2d1_2;
%     CSV_temp(:,4)=rsq2d1_10;
%    CSV_temp(:,4)=rsq2r1_2;
     CSV_temp(:,4)=rsq2r1_10;
    
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d1_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d2_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d10_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_df2_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_df10_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_rf2_new.csv');
    filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_rf10_new.csv');
    
    csvwrite(filename,CSV_temp);

    save('f20_1stloom_deltaNreg.mat');

%%

%%% should I try with the deltas or something else?

idx_temp2=idx_rsq_1stLoom_cleaned;
    
    CSV_temp=ROI_temp2.f20(idx_temp2,:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=rsquare_loom_short;
    
    %CSV_temp(:,5)=del_1_Nrsq;
    %CSV_temp(:,5)=del_2_Nrsq;
    CSV_temp(:,5)=rsq2d1_10;
    
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d1.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d2.csv');
    filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_df.csv');
    
    csvwrite(filename,CSV_temp);
    
    
   