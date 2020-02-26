%%% this script is to look at the max response to the looms of each dataset
%%% but filtering the ROIs to a first loom regression from
%%% testing_deltaNcorr.m script. 


%%% i need to add loading the idx_rsq and maybe some other things. 



load('All_More_BrainReg.mat','PerBrainRegions');

load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx', 'S_trim');


%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain'};

RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};

%% for f20

load('final_F20_step1.mat','idx_Fish_f20','ZS_f20');

%ZS_f20=ZS_f20.ZS_CN;

ZS_delta_per_fishNbrain_f20=struct;
fish=unique(idx_Fish_f20);
edges=[0:0.25:20];

for brain=1:length(RegionList)

    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_delta=intersect(temp_idx_fish,idx_rsq_1stLoom_cleaned);
temp_idx=intersect(temp_idx_delta,PerBrainRegions.f20.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_max_resp={NaN}; 
else
    
        for k=1:30

        temp_max_resp{1,k}= (max(ZS_f20(temp_idx,loom_moments{1,k})'))-(ZS_f20(temp_idx,Loomf20_onset_idx(k)+1))'; 
        
        [count_temp_max_resp{1,k},~]=histcounts(temp_max_resp{1,k},edges);
        end


end

ZS_delta_per_fishNbrain_f20.maxresponsePerloom.(RegionList{brain}){tempfish,1}=temp_max_resp;
ZS_delta_per_fishNbrain_f20.countMaxresponsePerloom.(RegionList{brain}){tempfish,1}=count_temp_max_resp;


end
end



%%

for brain=1:length(RegionList)
   
 for k=1:30
     
     temp_count=zeros(length(fish),length(edges)-1);
    for tempfish=1:length(fish)

    temp_count(tempfish,:)=ZS_delta_per_fishNbrain_f20.countMaxresponsePerloom.(RegionList{brain}){tempfish,1}{1,k};
    
    end
    
    temp_mean=mean(temp_count);
    ZS_delta_per_fishNbrain_f20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}){1,k}=temp_mean; 

end
end


%%% to check it worked
 figure('Position',[100 0 900 900]);
counter=1;
for brain=1:length(RegionList)

    
for i=[1 2 3 11]
    
    subplot(9,4,counter);
   plot(ZS_delta_per_fishNbrain_f20.MeancountMaxresponsePerloomPerfish.(RegionList{brain}){1,i}); title(RegionList{brain});
   
   counter=counter+1; 
end


sgtitle(RegionList{brain});

end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS ZS_delta_per_fishNbrain_f20


%% for f60 I cant until I also do a linear regression to the first loom. 
%%%for f60

load('final_F60_step1_2.mat','ZS_f60','idx_Fish_f60','ZS_short_F60');
%f60_cleaned_idxs=load('f60_cleaned_idxs.mat');


fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...

Stimuli=zeros(1,100);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=60;
for i=1%:10
    Stimuli(i,(idxStart+(i-1)*60):(idxStart+(i-1)*60)+size(GCaMP6,1)-1)=GCaMP6;
end


ModelResults=[];
parfor i=1:size(ZS_f60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(Stimuli',ZS_f60(i,1:100));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
     
end

rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them

figure; histogram(rsquare_loom);
std_rsq=std(rsquare_loom);
  

%%% what if I do 0.5?
idx_rsq=find(rsquare_loom>0.5 & rsquare_loom<1); %%%
figure; 
imagesc(ZS_f60(idx_rsq,:), [-0.5 4]);colormap hot



