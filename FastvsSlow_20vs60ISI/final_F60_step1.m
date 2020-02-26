
%%%% f60 final_step1

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

%%% for f60



F60data_CN=load('f60_postKmeans_CN_2.mat','ZS_CN','MatFiles');
ZS_f60=F60data_CN.('ZS_CN');
MatFiles_f60=F60data_CN.('MatFiles');



F60data=load('f60_r2050_CL4.mat','idx_Plane','idx_Fish');

idx_Plane_f60=F60data.('idx_Plane');
idx_Fish_f60=F60data.('idx_Fish');


%%% i am taking fish 37 away cause it itself has half the neurons (>10K) of the slope cluster.
ZS_f60(idx_Fish_f60==37,:)=[];

idx_Fish_allf60=idx_Fish_f60;

idx_Fish_f60(idx_Fish_allf60==37,:)=[];

idx_Plane_f60(idx_Fish_allf60==37,:)=[];


%%
%%% i am getting some weird results on the inhibition cluster using this
%%% data.. I will try to collect it directly. 

MatFiles=dir('*f60*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium

Noise=load(name, 'Noise');
Noise=Noise.Noise;

Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:);  %%%to combine the Calcium and Fitness variables (need to ask Gilles what Fitness is). Fitness here is the variable were we take the good calcium responses from the HPC analysis and pairthem with their index number.

GoodNoise=Noise(Fitness,:);



MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
% S=load(name, 'Spikes');
% S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;%%%because indexing in python is from 0 and matlab is at 1
%D=load(name, 'dFonF');
%D=D.dFonF;
GC=C(F,:);
%GS=S(F,:);
%GD=D(F,:);

    Noise=vertcat(Noise,N);
    GN=N(F,:);
    %Calcium=vertcat(Calcium,C);
    %DF=vertcat(DF,D);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC); %The fish 20+ are longer
    GoodNoise=vertcat(GoodNoise,GN);
    %GoodDF=vertcat(GoodDF,GD);
    %GoodSpikes=vertcat(GoodSpikes,GS);

%MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;%%%to get rid of vairables we will not use anymore

%%

GoodCalNoise=zeros(size(GoodNoise));
GoodCalNoise(:,:)=GoodCalcium+GoodNoise;

clear   GoodCalcium GoodNoise 

  save('f60_CN_GoodCalNoise_2.mat','-v7.3');


ZS_f60=zscore(GoodCalNoise,1,2); %%%to normalize the data
ZS_f60=detrend(ZS_f60')';%%% to Remove a linear trend from ZS (why???)




 clear Calcium Fitness GoodCalcium Noise GoodCalNoise
 
 save('f60_CN_2.mat','-v7.3');
 
 %%%the following is to create 1 column variables with the number of the
%%%fish and the slice number
Numbers=[1 [MatFiles.GoodNumber]]; %%%to take make a vector with the GoodNumber field
counter=1;
idx_Plane=nan(length(ZS_f60),1);%%%% to make an empty (with nans) one column variable the size of ZS
idx_Fish=nan(length(ZS_f60),1);%%%% to make an empty (with nans) one column variable the size of ZS
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %%%to get the number of the plane    
    idx_Plane(Numbers(i):Numbers(i+1))=Plane; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Plane   
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1}); %%%to get the number of the fish 
    idx_Fish(Numbers(i):Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore

idx_Plane_f60=idx_Plane;
idx_Fish_f60=idx_Fish; 

%%

rawregressF20=load('rawregressF20.mat','rawregress');
rawregressF20 = rawregressF20.('rawregress');

rawregressF20_CN=load('f20_CN_r2050_CL5_extra.mat','rawregress_CN');
rawregressF20_CN = rawregressF20_CN.('rawregress_CN');

rawregressF60=load('f60_r2050_CL4.mat','rawregress');
rawregressF60 = rawregressF60.('rawregress');



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(rawregressF20);

for i=1:length(rawregressF20)%% +2 %%%for the multisensregress
   
    subplot(4,2,counter);plot(rawregressF20(i,:));
     counter=counter+1;
end

rawregressF20(7,:)=rawregressF20_CN(1,:);


%%% to take the timepoints of f60 to fit it in f20
 ZS_short_F60=zeros(size(rawregressF60(1,:)));
 startpoints=[0,584,1168]; %%in seconds

 loom_times=[0,96,150,210,264,330,390,450,504,570]; %%in seconds
 loom_length=[52,18,20,18,22,20,20,18,22,14]; %%% in seconds

for p=1:3
    for k=1:10
    ZS_short_F60(startpoints(p)*2+loom_times(k)*2+1:startpoints(p)*2+loom_times(k)*2+loom_length(k)*2)=1;
    end
end


ZS_short_F60=find(ZS_short_F60==1);

%figure; plot(rawregressF60(1,ZS_short_F60)); %%% to check

%%

ModelResults_shortF60_all={};
for j=1:length(rawregressF20)
    
    
ModelResults_shortF60=[];
parfor i=1:size(ZS_f60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressF20(j,:)',ZS_f60(i,ZS_short_F60));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortF60(i).coef=mdl.Coefficients;
    ModelResults_shortF60(i).MSE=mdl.MSE;
    ModelResults_shortF60(i).Fitted=mdl.Fitted;
    ModelResults_shortF60(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortF60_all{j}=ModelResults_shortF60;

%rsquare_loom_shortF60=[ModelResults_shortF60.rsquared];


% figure;histogram(rsquare_loom_shortF60,'Normalization','probability');


%idx_rsq_test_f60short=find(rsquare_loom_shortF60>0.3 & rsquare_loom_shortF60<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
%proportion_f60=length(idx_rsq_test_f60short)/length(rsquare_loom_shortF60);

% figure;imagesc(ZS_f60(idx_rsq_test_f60short,:), [-0.5 4]);colormap hot

end




%%

idx_rsq_test_f60short_all={};
for j=1:size(rawregressF20)
    temp_rsq=[ModelResults_shortF60_all{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        %%% i tested 0.2, 0.25, 0.3 and 0.5. I think the best are 0.25 and 0.3 so far... 0.5 give nice clean
        %%% results but very few ROIs (12000), i think 0.25 is the limit of getting as much ROIs without gettig crap.
        %%% the shape and the distribution in planes and fish is very
        %%% similiar or equal between 0.25 and 0.3. but i get more ROIs
        %%% with 0.25.
        
        %%% with denoised data and 0.5 and the regressors together i get over 60K ROIs, which
        %%% is similar to what i get with 0.3 in the noised data and all
        %%% the regressors together. but by separating the regressors in a
        %%% loop i get better results but less ROIs
        
        %%% with 0.25, noised data and separated regressors i get 45K.
        %%% with 0.3 i get 35K and 0.2 is 55K.and just 12K with 0.5
        %%% I will check with denoised data
        %%% 0.5 gives me 35K
        
        %%% so it seems that 0.3 with noised data is  more or less the
        %%% equivalent to 0.5 with the denoised data.
    
idx_rsq_test_f60short_all{j}=idx_rsq_temp;

end

idx_rsq_test_f60short=horzcat(idx_rsq_test_f60short_all{:});
idx_rsq_test_f60short=unique(idx_rsq_test_f60short);


figure;imagesc(ZS_f60(idx_rsq_test_f60short,ZS_short_F60), [-0.5 4]);colormap hot


%%

%%% to clasify them with a correlation

Correlation_group_f60={};
counter=1;
for i=1:size(rawregressF20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_f60(idx_rsq_test_f60short),2)
        temp_corr=corrcoef(rawregressF20(i,:),ZS_f60(idx_rsq_test_f60short(idx),ZS_short_F60));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_f60{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_f60=[];
for n=1:size(rawregressF20,1)
Correlation_group_mat_f60(n,:)=Correlation_group_f60{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_f60=zeros(length(Correlation_group_mat_f60),1);
for i=1:length(Correlation_group_mat_f60)
    [~,I]=max(Correlation_group_mat_f60(:,i));
    High_corr_Nb_f60(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressF20,1);counter=1;
for i=1:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_F60_CN_each_2','-dpdf','-bestfit');

%%% i keep on getting weird ROIs in the inhibitory cluster. it seems that
%%% is from fish 3. They look like they are negative values due to
%%% movements....I confirmed that they belong to movements in the first slice of fish 3. I will need to clean them. 
%%% also, it seems that fish 47 has a loot of the ROIs in most of the
%%% clusters. specially broad fast hab (cluster 5 for f60)

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=1;counter=1;
for i=7%:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end


idx_temp=idx_rsq_test_f60short(find(High_corr_Nb_f60==7));
idx_fish3_inh=intersect(find(idx_Fish_f60==3),idx_temp);
%figure;imagesc(ZS_f60(idx_fish3_inhb,ZS_short_F60),[0 3]);
idx_mov_fish3_plane1=intersect(idx_fish3_inh,find(idx_Plane==1));
%figure;imagesc(ZS_f60(idx_mov_fish3_plane1,ZS_short_F60),[0 3]);
idx_todelete=find(ismember(idx_rsq_test_f60short,idx_mov_fish3_plane1));



idx_rsq_test_f60short2=idx_rsq_test_f60short;
High_corr_Nb_f60_2=High_corr_Nb_f60;


idx_rsq_test_f60short2(idx_todelete)=[];
High_corr_Nb_f60_2(idx_todelete)=[];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=1;counter=1;
for i=7%:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60_2==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short2(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short2(idx_temp)));%%% for the fish location
    counter=counter+4;
end




Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressF20,1);counter=1;
for i=1:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60_2==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short2(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short2(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_F60_CN_each_2_corrected','-dpdf','-bestfit');

%%

High_corr_Nb_f60_short=High_corr_Nb_f60_2;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_2==4))=2;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_2==5))=2;

gooodmaps=[1 2 6 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_f60_short==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short2(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short2(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_F60_CN_each_2','-dpdf','-bestfit');

%%

 save('final_F60_step1_2_short_correction.mat','idx_rsq_test_f60short2','High_corr_Nb_f60_2','idx_todelete','High_corr_Nb_f60_short','-v7.3');


%%% I also need to take away fish 47 as it has more than half of the
%%% responses for the fast hab broad cluster.
%%% i will do it through indexing to not repeat the whole analysis again.


load('final_F60_step1_2_short_correction.mat')

load('final_F60_step1_2.mat','MatFiles','ZS_f60','idx_Fish_f60','idx_Plane_f60','rawregressF20','ZS_short_F60');



%%% maybe not necesary
% ZS_f60(idx_Fish_f60==47,:)=[];
% 
% idx_Fish_allf60=idx_Fish_f60;
% 
% idx_Fish_f60(idx_Fish_allf60==47,:)=[];
% 
% idx_Plane_f60(idx_Fish_allf60==47,:)=[];




% %idx_temp=idx_rsq_test_f60short(find(High_corr_Nb_f60==7));
% idx_fish47_rsq=intersect(find(idx_Fish_allf60==47),idx_rsq_test_f60short2);
% %figure;imagesc(ZS_f60(idx_fish47_rsq,ZS_short_F60),[0 3]);
% %idx_mov_fish3_plane1=intersect(idx_fish47_rsq,find(idx_Plane==1));
% %figure;imagesc(ZS_f60(idx_mov_fish3_plane1,ZS_short_F60),[0 3]);
% idx_todelete=find(ismember(idx_fish47_rsq,idx_rsq_test_f60short2));

idx_todelete=find((idx_Fish_allf60(idx_rsq_test_f60short2)==47)>0);
idx_rsq_test_f60short3=idx_rsq_test_f60short2;
idx_rsq_test_f60short3(idx_todelete)=[];


High_corr_Nb_f60_3=High_corr_Nb_f60_2;

%idx_rsq_test_f60short3(idx_todelete)=[];
High_corr_Nb_f60_3(idx_todelete)=[];




Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressF20,1);counter=1;
for i=1:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60_3==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short3(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short3(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_F60_CN_each_3_corrected','-dpdf','-bestfit');

High_corr_Nb_f60_short=High_corr_Nb_f60_3;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_3==4))=2;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_3==5))=2;

gooodmaps=[1 2 6 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_f60_short==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short3(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short3(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_F60_CN_each_3','-dpdf','-bestfit');

 save('final_F60_step1_3_short_correction.mat','idx_rsq_test_f60short3','High_corr_Nb_f60_3','idx_todelete','High_corr_Nb_f60_short','-v7.3');