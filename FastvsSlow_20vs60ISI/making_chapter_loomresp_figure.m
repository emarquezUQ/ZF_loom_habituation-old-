%%%% this script is to get the ROIs, the mean and a rasterplot to make a
%%%% figure for the book chapter with loom respondent ROIs with some
%%%% emphasis in the pallium and Tegmentum. I will try with S20 cause i
%%%% didnt took any fish out. 

load('s20_cleaned_idxs.mat'); %%%% loading all the idx
load('final_S20_step1.mat','ZS_s20'); %%% loading the ZS


%%% here I am making a regressor for just the first loom. 
%%% based on the response 

figure;plot(mean(ZS_s20(idx_rsq_test_s20short_cleaned,:)));%%% so the start of the loom is at 60 and the spike falls at 100
figure;plot(mean(ZS_s20(clust_s20_CL4_cleaned.clust_s20_CL4_3_cleaned,:)));


Stimuli=zeros(1,100);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=60;
for i=1%:10
    Stimuli(i,(idxStart+(i-1)*60):(idxStart+(i-1)*60)+size(GCaMP6,1)-1)=GCaMP6;
end


%%
%%%% now I will do a linear regression but just to the ROIs that already
%%%% passed the first 0.3 original filter. 



    
ModelResults_shortS20=[];
parfor i=1:size(ZS_s20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(Stimuli(1,:)',ZS_s20(i,1:100));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS20(i).coef=mdl.Coefficients;
    ModelResults_shortS20(i).MSE=mdl.MSE;
    ModelResults_shortS20(i).Fitted=mdl.Fitted;
    ModelResults_shortS20(i).rsquared=mdl.Rsquared.Adjusted;
     
end

rsquare_loom=[ModelResults_shortS20.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them
idx_rsq=find(rsquare_loom>0.5 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; 
imagesc(ZS_s20(idx_rsq,:), [-0.5 4]);colormap hot %%%to plot them in a raster plot

idx_rsq2=find(rsquare_loom>0.3 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; 
imagesc(ZS_s20(idx_rsq2,:), [-0.5 4]);colormap hot %%%to plot them in a raster plot

idx_rsq3=find(rsquare_loom>0.4 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; 
imagesc(ZS_s20(idx_rsq3,:), [-0.5 4]);colormap hot %%%to plot them in a raster plot


figure;imagesc(ZS_s20(idx_rsq_test_s20short,:), [-0.5 4]);colormap hot %%% to compare



%%
%%%% to make the figures

figure;imagesc(ZS_s20(idx_rsq2,1:100),[min(min(mean(ZS_s20(idx_rsq2,1:100))))-0.5 max(max(mean(ZS_s20(idx_rsq2,1:100))))]);colormap hot
saveas(gcf,'rasterplot_s20_oneloom_forchapter2.emf');

figure;plot(mean(ZS_s20(idx_rsq2,1:100)));
saveas(gcf,'mean_s20_oneloom_forchapter2.emf');

figure;plot(Stimuli(1,:));
saveas(gcf,'Stimuli_s20_oneloom_forchapter.emf');


%%

%%% now to do the figure with the ROIs

%%% I need to get the Zbrain_Masks

load('Zbrain_Masks.mat');

load('All_More_BrainReg2.mat');


%%% this is to make a matrix with the location of the main structures 
Zbrain_brainMask=vertcat(Zbrain_Masks{[76 113 259 274 294],3}); %%% 78 is the eyes
 
Zbrain_brainMask=unique(Zbrain_brainMask,'rows');

figure;scatter(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),'.');
%figure;scatter3(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),Zbrain_brainMask(:,3),'.');


figure;scatter(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),'filled','k','MarkerFaceAlpha',.01);
hold on;
scatter(ROI_temp2.s20(idx_rsq2,1),ROI_temp2.s20(idx_rsq2,2),2,'filled','w');
hold on;
scatter(ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),2),5,'filled','r');
hold on;
scatter(ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq2),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq2),2),5,'filled','g');

saveas(gcf,'BrainNrois_s20_oneloom_forchapter.emf');


%% 

%%% if i use the fast hab ones. 

figure;scatter(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),'filled','k','MarkerFaceAlpha',.01);
hold on;
scatter(ROI_temp2.s20(clust_s20_CL4_cleaned.clust_s20_CL4_1_cleaned,1),ROI_temp2.s20(clust_s20_CL4_cleaned.clust_s20_CL4_1_cleaned,2),2,'filled','w');
hold on;
scatter(ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,clust_s20_CL4_cleaned.clust_s20_CL4_1_cleaned),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,clust_s20_CL4_cleaned.clust_s20_CL4_1_cleaned),2),5,'filled','r');
hold on;
scatter(ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,clust_s20_CL4_cleaned.clust_s20_CL4_1_cleaned),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,clust_s20_CL4_cleaned.clust_s20_CL4_1_cleaned),2),5,'filled','g');

saveas(gcf,'BrainNrois_s20_oneloom_forchapter_CL4.emf');

%%

%%% for a 3d image

%%%% with R2=0.3

patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(ROI_temp2.s20(idx_rsq2,1),ROI_temp2.s20(idx_rsq2,2),ROI_temp2.s20(idx_rsq2,3),2,'filled','k');
hold on;
scatter3(ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),2),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),3),5,'filled','r');
hold on;
scatter3(ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq2),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq2),2),ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq2),3),5,'filled','g');

saveas(gcf,'BrainNrois_s20_oneloom_forchapter3D.emf');

%%% now just one of the pallium in frontal view

%%% to test if it works

patch(BrainRegions3D.Pallium.poly,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),2),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq2),3),5,'o','filled','r');
%set(gca,'Color','k'); %%% color of background to black

saveas(gcf,'BrainNrois_s20_oneloom_forchapter3D_pallium.emf');



%%

%%%% with R2=0.4

%%% for a 3d image

patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(ROI_temp2.s20(idx_rsq3,1),ROI_temp2.s20(idx_rsq3,2),ROI_temp2.s20(idx_rsq3,3),2,'filled','k');
hold on;
scatter3(ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq3),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq3),2),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq3),3),5,'filled','r');
hold on;
scatter3(ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq3),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq3),2),ROI_temp2.s20(intersect(PerBrainRegions.s20.Tegmentum.idx,idx_rsq3),3),5,'filled','g');

%saveas(gcf,'BrainNrois_s20_oneloom_forchapter3D_R040.emf');

%%% now just one of the pallium in frontal view

%%% to test if it works

patch(BrainRegions3D.Pallium.poly,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq3),1),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq3),2),ROI_temp2.s20(intersect(PerBrainRegions.s20.Pallium.idx,idx_rsq3),3),10,'o','filled','r');
%set(gca,'Color','k'); %%% color of background to black

%saveas(gcf,'BrainNrois_s20_oneloom_forchapter3D_pallium_R040.emf');