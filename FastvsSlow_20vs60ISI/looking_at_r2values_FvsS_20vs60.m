%%%% this script is to look at the r2 values of the 4 datasets


%%% it seems that the s60 has a weird result... far less neurons than the
%%% others. 


cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab


rsquare_loom_s60=load('s60_r2050_CL3.mat','rsquare_loom');
rsquare_loom_s60 = rsquare_loom_s60.('rsquare_loom');


figure;histogram(rsquare_loom_s60,'Normalization','probability');


idx_rsq_test_s60=find(rsquare_loom_s60>0.5 & rsquare_loom_s60<1);
length(idx_rsq_test_s60)/length(rsquare_loom_s60)


rsquare_loom_f60=load('f60_r2050_CL4.mat','rsquare_loom');
rsquare_loom_f60 = rsquare_loom_f60.('rsquare_loom');


figure;histogram(rsquare_loom_f60,'Normalization','probability');


idx_rsq_test_f60=find(rsquare_loom_f60>0.5 & rsquare_loom_f60<1);
length(idx_rsq_test_f60)/length(rsquare_loom_f60)



rsquare_loom_s20=load('s20_r2050_CL6.mat','rsquare_loom');
rsquare_loom_s20 = rsquare_loom_s20.('rsquare_loom');


figure;histogram(rsquare_loom_s20,'Normalization','probability');

idx_rsq_test_s20=find(rsquare_loom_s20>0.5 & rsquare_loom_s20<1);
length(idx_rsq_test_s20)/length(rsquare_loom_s20)


rsquare_loom_f20=load('f20_r2050_CL5.mat','rsquare_loom');
rsquare_loom_f20 = rsquare_loom_f20.('rsquare_loom');


figure;histogram(rsquare_loom_f20,'Normalization','probability');

idx_rsq_test_f20=find(rsquare_loom_f20>0.5 & rsquare_loom_f20<1);
length(idx_rsq_test_f20)/length(rsquare_loom_f20)


%%


%%% here i will try to use the regressors from the s20 to see if I can find
%%% more ROIs in s60


%%% NOTE: it worked!!! i got 2 times more ROIs that passed the 0.5
%%% threshold of the linear regression. Probably i was underclustering. 
%%% Gilles suggests that I should do this for f20-f60 too. 

S60data=load('s60_r2050_CL3.mat','ZS','idx_rsq','idx_coef_rsq','idxKmeans1_coef_rsq','idxKmeans_ZS','idx_Plane','idx_Fish','goodmaps');
ZS_s60=S60data.('ZS');
idx_rsq_s60=S60data.('idx_rsq');
idx_coef_rsq_s60=S60data.('idx_coef_rsq');
idxKmeans1_coef_rsq_s60=S60data.('idxKmeans1_coef_rsq');
idxKmeans_ZS_s60=S60data.('idxKmeans_ZS');
idx_Plane_s60=S60data.('idx_Plane');
idx_Fish_s60=S60data.('idx_Fish');
goodmaps_s60=S60data.('goodmaps');



S20data=load('s20_r2050_CL6.mat','ZS','idx_rsq','idx_coef_rsq','idxKmeans1_coef_rsq','idxKmeans_ZS','idx_Plane','idx_Fish','goodmaps');
ZS_s20=S20data.('ZS');
idx_rsq_s20=S20data.('idx_rsq');
idx_coef_rsq_s20=S20data.('idx_coef_rsq');
idxKmeans1_coef_rsq_s20=S20data.('idxKmeans1_coef_rsq');
idxKmeans_ZS_s20=S20data.('idxKmeans_ZS');
idx_Plane_s20=S20data.('idx_Plane');
idx_Fish_s20=S20data.('idx_Fish');
goodmaps_s20=S20data.('goodmaps');


rawregressS20=load('rawregressS20.mat','rawregress');
rawregressS20 = rawregressS20.('rawregress');

rawregressS60=load('rawregressS60.mat','rawregress');
rawregressS60 = rawregressS60.('rawregress');

%%% to take the timepoints of s60 to fit it in s20
 ZS_short_S60=zeros(size(rawregressS60(1,:)));
 startpoints=[0,586,1172]; %%in seconds

 loom_times=[0,96,150,210,264,330,390,450,504,570]; %%in seconds
 loom_length=[52,18,20,18,22,20,20,18,22,16]; %%% in seconds

for p=1:3
    for k=1:10
    ZS_short_S60(startpoints(p)*2+loom_times(k)*2+1:startpoints(p)*2+loom_times(k)*2+loom_length(k)*2)=1;
    end
end


ZS_short_S60=find(ZS_short_S60==1);

% figure; plot(rawregressS60(2,:));
% 
% figure; plot(rawregressS60(2,ZS_short_S60));
% 
% length(ZS_short_S60)
% 
% length(rawregressS20)



%figure;plot(ZS_s60(idx_rsq_s60(50),ZS_short_S60)); %%% to check


ModelResults_shortS60=[];
parfor i=1:size(ZS_s60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressS20',ZS_s60(i,ZS_short_S60));
    %mdl=stepwiselm(rawregress',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS60(i).coef=mdl.Coefficients;
    ModelResults_shortS60(i).MSE=mdl.MSE;
    ModelResults_shortS60(i).Fitted=mdl.Fitted;
    ModelResults_shortS60(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom_shortS60=[ModelResults_shortS60.rsquared];


figure;histogram(rsquare_loom_shortS60,'Normalization','probability');


idx_rsq_test_s60short=find(rsquare_loom_shortS60>0.5 & rsquare_loom_shortS60<1);
length(idx_rsq_test_s60short)/length(rsquare_loom_shortS60)

figure;imagesc(ZS_s60(idx_rsq_test_s60short,:), [-0.5 4]);colormap hot



figure;imagesc(ZS_s60(idx_rsq_s60,:), [-0.5 4]);colormap hot



%%%this is to at the clusters that i had before on the s60.

%%% i am getting less ROIs in this graph... It is because i only used 4
%%% of the clusters that got out of the first Kmeans (some of them were
%%% not in all fish). 

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(unique(goodmaps_s60));counter=1;
for i=unique(goodmaps_s60)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq_s60==i);
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_coef_rsq_s60(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_coef_rsq_s60(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_coef_rsq_s60(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_coef_rsq_s60(idx_temp)));%%% for the fish location
    counter=counter+4;
end


%%%this is look to at the clusters that i had before on the s20.
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(unique(goodmaps_s20));counter=1;
for i=unique(goodmaps_s20)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_coef_rsq_s20(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_coef_rsq_s20(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_coef_rsq_s20(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_coef_rsq_s20(idx_temp)));%%% for the fish location
    counter=counter+4;
end


%%%

%%% checking at the slop hab recovery of the old ones. i will normilize
%%% them too.

%%% comparing the means now and before of the slope habcluster... it does seem to be a small
%%% diference in the recovery but normalzing the way i am doing it doesnt
%%% seem to work very well. 
idx_temp1=find(idxKmeans1_coef_rsq_s60==8);
idx_temp2=find(idxKmeans1_coef_rsq_s20==41);



figure;plot(mean(ZS_s60(idx_coef_rsq_s60(idx_temp1),ZS_short_S60),1));
hold on; plot(mean(ZS_s20(idx_coef_rsq_s20(idx_temp2),:),1));

hold off;

norm_mean_slopehab_s60=mean(ZS_s60(idx_coef_rsq_s60(idx_temp1),ZS_short_S60),1);

for i=1:length(norm_mean_slopehab_s60)
    norm_mean_slopehab_s60(1,:)=norm_mean_slopehab_s60(1,:)-mean(norm_mean_slopehab_s60(1,1:5));
    norm_mean_slopehab_s60(1,:)=norm_mean_slopehab_s60(1,:)/max(norm_mean_slopehab_s60(1,:));

end

%figure;plot(norm_mean_slopehab_s60);


norm_mean_slopehab_s20=mean(ZS_s20(idx_coef_rsq_s20(idx_temp2),:),1);

for i=1:length(norm_mean_slopehab_s20)
    norm_mean_slopehab_s20(1,:)=norm_mean_slopehab_s20(1,:)-mean(norm_mean_slopehab_s20(1,1:5));
    norm_mean_slopehab_s20(1,:)=norm_mean_slopehab_s20(1,:)/max(norm_mean_slopehab_s20(1,:));

end

%figure;plot(norm_mean_slopehab_s20);

figure;plot(norm_mean_slopehab_s60);
hold on; plot(norm_mean_slopehab_s20);

hold off;

figure;
subplot(2,1,1); plot(norm_mean_slopehab_s60);
subplot(2,1,2); plot(norm_mean_slopehab_s20);



%%

%%% to clasify them with a correlation

Correlation_group_s60={};
counter=1;
for i=1:size(rawregressS20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s60(idx_rsq_test_s60short),2)
        temp_corr=corrcoef(rawregressS20(i,:),ZS_s60(idx_rsq_test_s60short(idx),ZS_short_S60));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_s60{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_s60=[];
for n=1:size(rawregressS20,1)
Correlation_group_mat_s60(n,:)=Correlation_group_s60{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_s60=zeros(length(Correlation_group_mat_s60),1);
for i=1:length(Correlation_group_mat_s60)
    [~,I]=max(Correlation_group_mat_s60(:,i));
    High_corr_Nb_s60(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(High_corr_Nb_s60==i);
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_rsq_test_s60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_rsq_test_s60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%%


%%% now with the s20 again


ModelResults_newS20=[];
parfor i=1:size(ZS_s20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressS20',ZS_s20(i,:));
    %mdl=stepwiselm(rawregress',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_newS20(i).coef=mdl.Coefficients;
    ModelResults_newS20(i).MSE=mdl.MSE;
    ModelResults_newS20(i).Fitted=mdl.Fitted;
    ModelResults_newS20(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom_newS20=[ModelResults_newS20.rsquared];
figure;histogram(rsquare_loom_newS20,'Normalization','probability');


figure;histogram(rsquare_loom_s20,'Normalization','probability');


idx_rsq_test_new_s20=find(rsquare_loom_newS20>0.5 & rsquare_loom_newS20<1);
length(idx_rsq_test_new_s20)/length(rsquare_loom_newS20)



%%

%%% to clasify them with a correlation

Correlation_group_s20={};
counter=1;
for i=1:size(rawregressS20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s20(idx_rsq_test_new_s20),2)
        temp_corr=corrcoef(rawregressS20(i,:),ZS_s20(idx_rsq_test_new_s20(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_s20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_s20=[];
for n=1:size(rawregressS20,1)
Correlation_group_mat_s20(n,:)=Correlation_group_s20{n}(n,:);

end

%%%what if I filter looking for the max correlation?
%%% it seems to be working!!!!

High_corr_Nb_s20=zeros(length(Correlation_group_mat_s20),1);
for i=1:length(Correlation_group_mat_s20)
    [~,I]=max(Correlation_group_mat_s20(:,i));
    High_corr_Nb_s20(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(High_corr_Nb_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_rsq_test_new_s20(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_rsq_test_new_s20(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_rsq_test_new_s20(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_rsq_test_new_s20(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%%

%%% lets try to put it together
%%% in looks very interesting. 


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp1=find(High_corr_Nb_s20==i);
    idx_temp2=find(High_corr_Nb_s60==i);
    subplot(rows,5,counter);plot(mean(ZS_s20(idx_rsq_test_new_s20(idx_temp1),:),1)); hold on; plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp2),ZS_short_S60),1)); hold off;%%%to plot the mean
    subplot(rows,5,counter+1);imagesc(ZS_s20(idx_rsq_test_new_s20(idx_temp1),:),[0 3]);%%%for the raster plot
    subplot(rows,5,counter+2);imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp2),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,5,counter+3);histogram(idx_Plane_s20(idx_rsq_test_new_s20(idx_temp1))); hold on; histogram(idx_Plane_s60(idx_rsq_test_s60short(idx_temp2))); hold off;%%%for the plane location   
    subplot(rows,5,counter+4);histogram(idx_Fish_s20(idx_rsq_test_new_s20(idx_temp1))); hold on; histogram(idx_Fish_s60(idx_rsq_test_s60short(idx_temp2))); hold off;%%% for the fish location
    counter=counter+5;
end


%print(gcf,'multigraph_s20ns60','-dpdf','-bestfit');

%%% this is just for the means. sadly it doesnt look as impresive as i
%%% thought it would. Maybe I would need to look at the variability.

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;%rows=length(rawregressS20);

for i=1:length(rawregressS20)
   idx_temp1=find(High_corr_Nb_s20==i);
    idx_temp2=find(High_corr_Nb_s60==i);
    subplot(3,2,counter);plot(mean(ZS_s20(idx_rsq_test_new_s20(idx_temp1),:),1)); hold on; plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp2),ZS_short_S60),1)); hold off;%%%to plot the mean
    
     counter=counter+1;
end


%%% comparing the means now and before of the slope habcluster... it doest seem to be a small
%%% diference in the recovery but is not super impresive




idx_temp1=find(High_corr_Nb_s60==5);
idx_temp2=find(High_corr_Nb_s20==5);     

figure; plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp1),ZS_short_S60),1));     
hold on; plot(mean(ZS_s20(idx_rsq_test_new_s20(idx_temp2),:),1));
hold off;



%%% if i normalize it

norm_mean_slopehab_s60_new=mean(ZS_s60(idx_rsq_test_s60short(idx_temp1),ZS_short_S60),1);

for i=1:length(norm_mean_slopehab_s60_new)
    norm_mean_slopehab_s60_new(1,:)=norm_mean_slopehab_s60_new(1,:)-mean(norm_mean_slopehab_s60_new(1,1:5));
    norm_mean_slopehab_s60_new(1,:)=norm_mean_slopehab_s60_new(1,:)/max(norm_mean_slopehab_s60_new(1,:));

end

%figure;plot(norm_mean_slopehab_s60);


norm_mean_slopehab_s20_new=mean(ZS_s20(idx_rsq_test_new_s20(idx_temp2),:),1);

for i=1:length(norm_mean_slopehab_s20_new)
    norm_mean_slopehab_s20_new(1,:)=norm_mean_slopehab_s20_new(1,:)-mean(norm_mean_slopehab_s20_new(1,1:5));
    norm_mean_slopehab_s20_new(1,:)=norm_mean_slopehab_s20_new(1,:)/max(norm_mean_slopehab_s20_new(1,:));

end

%figure;plot(norm_mean_slopehab_s20);


%%% there is a kind of problem. as I am using the s20 regressor with a
%%% higher recovery the new s60 mean gets a higher recovery than what it use to be. although it seems that the difference between ISIs is still there. 

figure;
subplot(2,1,1); plot(norm_mean_slopehab_s60);hold on;
subplot(2,1,1); plot(norm_mean_slopehab_s60_new); hold off;
subplot(2,1,2); plot(norm_mean_slopehab_s20);hold on;
subplot(2,1,2);plot(norm_mean_slopehab_s20_new);hold off;



%%
%%% Gilles suggests that I should do this for f20-f60 too. 

%%% it gives similar results than before. but the good thing now is that I
%%% will be able to look for the inhibitors and sound responding ROIs


%%% for f60

F60data=load('f60_r2050_CL4.mat','ZS','idx_rsq','idx_coef_rsq','idxKmeans1_coef_rsq','idxKmeans_ZS','idx_Plane','idx_Fish','goodmaps');
ZS_f60=F60data.('ZS');
idx_rsq_f60=F60data.('idx_rsq');
idx_coef_rsq_f60=F60data.('idx_coef_rsq');
idxKmeans1_coef_rsq_f60=F60data.('idxKmeans1_coef_rsq');
idxKmeans_ZS_f60=F60data.('idxKmeans_ZS');
idx_Plane_f60=F60data.('idx_Plane');
idx_Fish_f60=F60data.('idx_Fish');
goodmaps_f60=F60data.('goodmaps');

%%% i am taking fish 37 away cause it itself has half the neurons (>10K) of the slope cluster.
ZS_f60(idx_Fish_f60==37,:)=[];

ZS_f60_old=F60data.('ZS');




%%% for f20

F20data=load('f20_r2050_CL5.mat','ZS','idx_rsq','idx_coef_rsq','idxKmeans1_coef_rsq','idxKmeans_ZS','idx_Plane','idx_Fish','goodmaps','GoodBetas_ZS2');
ZS_f20=F20data.('ZS');
idx_rsq_f20=F20data.('idx_rsq');
idx_coef_rsq_f20=F20data.('idx_coef_rsq');
idxKmeans1_coef_rsq_f20=F20data.('idxKmeans1_coef_rsq');
idxKmeans_ZS_f20=F20data.('idxKmeans_ZS');
idx_Plane_f20=F20data.('idx_Plane');
idx_Fish_f20=F20data.('idx_Fish');
goodmaps_f20=F20data.('goodmaps');
GoodBetas_ZS2_f20=F20data.('GoodBetas_ZS2');

rawregressF20=load('rawregressF20.mat','rawregress');
rawregressF20 = rawregressF20.('rawregress');

PostregressF20=load('f20_r2050_CL5.mat','Cmap_ZS_rsq');
PostregressF20 = PostregressF20.('Cmap_ZS_rsq');

rawregressF60=load('f60_r2050_CL4.mat','rawregress');
rawregressF60 = rawregressF60.('rawregress');



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(PostregressF20);

for i=1:length(PostregressF20)%% +2 %%%for the multisensregress
   
    subplot(4,2,counter);plot(PostregressF20(i,:));
     counter=counter+1;
end




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

%figure; plot(rawregressF60(1,:));

%figure; plot(rawregressF20(1,:));

%figure; plot(rawregressF60(1,ZS_short_F60));

%length(ZS_short_F60)

%length(rawregressF20)

%length(rawregressF60)

%figure;plot(ZS_f60(idx_rsq_f60(50),ZS_short_F60)); %%% to check


ModelResults_shortF60=[];
parfor i=1:size(ZS_f60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(PostregressF20',ZS_f60(i,ZS_short_F60));
    %mdl=stepwiselm(rawregress',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortF60(i).coef=mdl.Coefficients;
    ModelResults_shortF60(i).MSE=mdl.MSE;
    ModelResults_shortF60(i).Fitted=mdl.Fitted;
    ModelResults_shortF60(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom_shortF60=[ModelResults_shortF60.rsquared];
figure;histogram(rsquare_loom_shortF60,'Normalization','probability');

rsquare_loom_f20=load('f20_r2050_CL5.mat','rsquare_loom');
rsquare_loom_f20 = rsquare_loom_f20.('rsquare_loom');
figure;histogram(rsquare_loom_f20,'Normalization','probability');


idx_rsq_test_f60short=find(rsquare_loom_shortF60>0.5 & rsquare_loom_shortF60<1);
length(idx_rsq_test_f60short)/length(rsquare_loom_shortF60)

%figure;imagesc(ZS_f60(idx_coef_rsq_f60,ZS_short_F60), [-0.5 4]);colormap hot
%figure;imagesc(ZS_f60(idx_rsq_test_f60short,ZS_short_F60), [-0.5 4]);colormap hot




%%%this is to at the clusters that i had before on the f60.

%%% i am getting less ROIs in this graph... It is because i only used 4
%%% of the clusters that got out of the first Kmeans (some of them were
%%% not in all fish). 

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(unique(goodmaps_f60));counter=1;
for i=unique(goodmaps_f60)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq_f60==i);
    subplot(rows,4,counter);plot(mean(ZS_f60_old(idx_coef_rsq_f60(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60_old(idx_coef_rsq_f60(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_coef_rsq_f60(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_coef_rsq_f60(idx_temp)));%%% for the fish location
    counter=counter+4;
end


%%%this is look to at the clusters that i had before on the f20.
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(unique(GoodBetas_ZS2_f20));counter=1;
for i=unique(GoodBetas_ZS2_f20)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq_f20==i);
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_coef_rsq_f20(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_coef_rsq_f20(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_coef_rsq_f20(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_coef_rsq_f20(idx_temp)));%%% for the fish location
    counter=counter+4;
end


%%% checking at the slop hab recovery of the old ones. i will normilize
%%% them too.

%%% comparing the means now and before of the slope habcluster... it does seem to be a small
%%% diference in the recovery but normalzing the way i am doing it doesnt
%%% seem to work very well. 
idx_temp1=find(idxKmeans1_coef_rsq_f60==33);
idx_temp2=find(idxKmeans1_coef_rsq_f20==1);



figure;plot(mean(ZS_f60_old(idx_coef_rsq_f60(idx_temp1),ZS_short_F60),1));
hold on; plot(mean(ZS_f20(idx_coef_rsq_f20(idx_temp2),:),1));

hold off;

norm_mean_slopehab_f60=mean(ZS_f60_old(idx_coef_rsq_f60(idx_temp1),ZS_short_F60),1);

for i=1:length(norm_mean_slopehab_f60)
    norm_mean_slopehab_f60(1,:)=norm_mean_slopehab_f60(1,:)-mean(norm_mean_slopehab_f60(1,1:5));
    norm_mean_slopehab_f60(1,:)=norm_mean_slopehab_f60(1,:)/max(norm_mean_slopehab_f60(1,:));

end

%figure;plot(norm_mean_slopehab_f60);


norm_mean_slopehab_f20=mean(ZS_f20(idx_coef_rsq_f20(idx_temp2),:),1);

for i=1:length(norm_mean_slopehab_f20)
    norm_mean_slopehab_f20(1,:)=norm_mean_slopehab_f20(1,:)-mean(norm_mean_slopehab_f20(1,1:5));
    norm_mean_slopehab_f20(1,:)=norm_mean_slopehab_f20(1,:)/max(norm_mean_slopehab_f20(1,:));

end

%figure;plot(norm_mean_slopehab_f20);

figure;plot(norm_mean_slopehab_f60);
hold on; plot(norm_mean_slopehab_f20);

hold off;

figure;
subplot(2,1,1); plot(norm_mean_slopehab_f60);
subplot(2,1,2); plot(norm_mean_slopehab_f20);




%%

%%% to clasify them with a correlation

Correlation_group_f60={};
counter=1;
for i=1:size(PostregressF20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_f60(idx_rsq_test_f60short),2)
        temp_corr=corrcoef(PostregressF20(i,:),ZS_f60(idx_rsq_test_f60short(idx),ZS_short_F60));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_f60{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_f60=[];
for n=1:size(PostregressF20,1)
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
rows=size(PostregressF20,1);counter=1;
for i=1:size(PostregressF20,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(High_corr_Nb_f60==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end




%%

%%% now with the f20 again


ModelResults_newF20=[];
parfor i=1:size(ZS_f20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(PostregressF20',ZS_f20(i,:));
    %mdl=stepwiselm(rawregress',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_newF20(i).coef=mdl.Coefficients;
    ModelResults_newF20(i).MSE=mdl.MSE;
    ModelResults_newF20(i).Fitted=mdl.Fitted;
    ModelResults_newF20(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom_newF20=[ModelResults_newF20.rsquared];
figure;histogram(rsquare_loom_newF20,'Normalization','probability');


figure;histogram(rsquare_loom_f20,'Normalization','probability');


idx_rsq_test_new_f20=find(rsquare_loom_newF20>0.5 & rsquare_loom_newF20<1);
length(idx_rsq_test_new_f20)/length(rsquare_loom_newF20)



%%

%%% to clasify them with a correlation

Correlation_group_f20={};
counter=1;
for i=1:size(PostregressF20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_f20(idx_rsq_test_new_f20),2)
        temp_corr=corrcoef(PostregressF20(i,:),ZS_f20(idx_rsq_test_new_f20(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_f20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_f20=[];
for n=1:size(PostregressF20,1)
Correlation_group_mat_f20(n,:)=Correlation_group_f20{n}(n,:);

end

%%%what if I filter looking for the max correlation?
%%% it seems to be working!!!!

High_corr_Nb_f20=zeros(length(Correlation_group_mat_f20),1);
for i=1:length(Correlation_group_mat_f20)
    [~,I]=max(Correlation_group_mat_f20(:,i));
    High_corr_Nb_f20(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(PostregressF20,1);counter=1;
for i=1:size(PostregressF20,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(High_corr_Nb_f20==i);
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_rsq_test_new_f20(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_rsq_test_new_f20(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_rsq_test_new_f20(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_rsq_test_new_f20(idx_temp)));%%% for the fish location
    counter=counter+4;
end


%%

%%% lets try to put it together
%%% in looks very interesting. but i might need to take one of the f60 fish out (fish 38 or so) as its half the neurons for the slope hab cluster.


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(PostregressF20,1);counter=1;
for i=1:size(PostregressF20,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp1=find(High_corr_Nb_f20==i);
    idx_temp2=find(High_corr_Nb_f60==i);
    subplot(rows,5,counter);plot(mean(ZS_f20(idx_rsq_test_new_f20(idx_temp1),:),1)); hold on; plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp2),ZS_short_F60),1)); hold off;%%%to plot the mean
    subplot(rows,5,counter+1);imagesc(ZS_f20(idx_rsq_test_new_f20(idx_temp1),:),[0 3]);%%%for the raster plot
    subplot(rows,5,counter+2);imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp2),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,5,counter+3);histogram(idx_Plane_f20(idx_rsq_test_new_f20(idx_temp1))); hold on; histogram(idx_Plane_f60(idx_rsq_test_f60short(idx_temp2))); hold off;%%%for the plane location   
    subplot(rows,5,counter+4);histogram(idx_Fish_f20(idx_rsq_test_new_f20(idx_temp1))); hold on; histogram(idx_Fish_f60(idx_rsq_test_f60short(idx_temp2))); hold off;%%% for the fish location
    counter=counter+5;
end


%%% this is just for the means. sadly it doesnt look as impresive as i
%%% thought it would. Maybe I would need to look at the variability.

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;%rows=length(PostregressF20);

for i=1:length(PostregressF20)
   idx_temp1=find(High_corr_Nb_f20==i);
    idx_temp2=find(High_corr_Nb_f60==i);
    subplot(3,2,counter);plot(mean(ZS_f20(idx_rsq_test_new_f20(idx_temp1),:),1)); hold on; plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp2),ZS_short_F60),1)); hold off;%%%to plot the mean
    
     counter=counter+1;
end

%%% comparing the means now and before of the slope habcluster... it doest seem to be a small
%%% diference in the recovery but is not very impresive
idx_temp1=find(High_corr_Nb_f20==1);
idx_temp2=find(High_corr_Nb_f60==1);
idx_temp3=find(idxKmeans1_coef_rsq_f60==33);
figure;plot(mean(ZS_f20(idx_rsq_test_new_f20(idx_temp1),:),1));
hold on; plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp2),ZS_short_F60),1)); 
hold on; plot(mean(ZS_f60_old(idx_coef_rsq_f60(idx_temp3),ZS_short_F60),1));
hold off;


idx_temp1=find(High_corr_Nb_f60==1);
idx_temp2=find(High_corr_Nb_f20==1);     

figure; plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp1),ZS_short_F60),1));     
hold on; plot(mean(ZS_f20(idx_rsq_test_new_f20(idx_temp2),:),1));
hold off;



%%% if i normalize it

norm_mean_slopehab_f60_new=mean(ZS_f60(idx_rsq_test_f60short(idx_temp1),ZS_short_F60),1);

for i=1:length(norm_mean_slopehab_f60_new)
    norm_mean_slopehab_f60_new(1,:)=norm_mean_slopehab_f60_new(1,:)-mean(norm_mean_slopehab_f60_new(1,1:5));
    norm_mean_slopehab_f60_new(1,:)=norm_mean_slopehab_f60_new(1,:)/max(norm_mean_slopehab_f60_new(1,:));

end

%figure;plot(norm_mean_slopehab_f60);


norm_mean_slopehab_f20_new=mean(ZS_f20(idx_rsq_test_new_f20(idx_temp2),:),1);

for i=1:length(norm_mean_slopehab_f20_new)
    norm_mean_slopehab_f20_new(1,:)=norm_mean_slopehab_f20_new(1,:)-mean(norm_mean_slopehab_f20_new(1,1:5));
    norm_mean_slopehab_f20_new(1,:)=norm_mean_slopehab_f20_new(1,:)/max(norm_mean_slopehab_f20_new(1,:));

end

%figure;plot(norm_mean_slopehab_f20);


%%% there is a kind of problem. as I am using the f20 regressor with a
%%% higher recovery the new f60 mean gets a higher recovery than what it use to be. although it seems that the difference between ISIs is still there. 

figure;
subplot(2,1,1); plot(norm_mean_slopehab_f60);hold on;
subplot(2,1,1); plot(norm_mean_slopehab_f60_new); hold off;
subplot(2,1,2); plot(norm_mean_slopehab_f20);hold on;
subplot(2,1,2);plot(norm_mean_slopehab_f20_new);hold off;


figure;
 plot(norm_mean_slopehab_f60);hold on;
 plot(norm_mean_slopehab_f60_new); hold on;
 plot(norm_mean_slopehab_f20);hold on;
 plot(norm_mean_slopehab_f20_new);hold on;


 %%
 %%% trying to find the way to quatify the recovery difference
 [psor,lsor]=findpeaks(norm_mean_slopehab_f60(375:575),'SortStr','descend','MinPeakDistance',40,'MinPeakHeight',0.05); 
 figure;findpeaks(norm_mean_slopehab_f60(375:575),'SortStr','descend','MinPeakDistance',40,'MinPeakHeight',0.1); text (lsor+.02,psor,num2str((1:numel(psor))'));
 
 
 %%% i will leave it for later
 
 
 