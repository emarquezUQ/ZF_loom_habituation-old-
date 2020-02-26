
%%%% s20 final_step1

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab


%%% now the ZS with the noise, matfiles,regressors, planes and fish. i think at the
%%% moment i am not excluding any of the Slow loom fish...
load('inhib_s20_regress_200CL.mat','inhib_s20_regress');


S20data_CN=load('s20_postKmeans_CN.mat','ZS_CN','MatFiles');
ZS_s20=S20data_CN.('ZS_CN');
MatFiles_s20=S20data_CN.('MatFiles');

S20data=load('s20_r2050_CL6.mat','idx_Plane','idx_Fish');

idx_Plane_s20=S20data.('idx_Plane');
idx_Fish_s20=S20data.('idx_Fish');


rawregressS20=load('rawregressS20.mat','rawregress');
rawregressS20 = rawregressS20.('rawregress');

rawregressS60=load('rawregressS60.mat','rawregress');
rawregressS60 = rawregressS60.('rawregress');


%%


%%% this is for a linear regression using the s20 denoised data regressors
%%% from the 50CL + the inhibition cluster

% counter=1;
% figure;
% for i=1:length(rawregressS20)
% subplot(4,2,counter);plot(rawregressS20(i,:));
% counter=counter+1;
% end


rawregressS20(7,:)=inhib_s20_regress;


ModelResults_shortS20_all={};
for j=1:length(rawregressS20)
    
    
ModelResults_shortS20=[];
parfor i=1:size(ZS_s20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressS20(j,:)',ZS_s20(i,:));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS20(i).coef=mdl.Coefficients;
    ModelResults_shortS20(i).MSE=mdl.MSE;
    ModelResults_shortS20(i).Fitted=mdl.Fitted;
    ModelResults_shortS20(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortS20_all{j}=ModelResults_shortS20;

%rsquare_loom_shortS20=[ModelResults_shortS20.rsquared];


% figure;histogram(rsquare_loom_shortS20,'Normalization','probability');


%idx_rsq_test_s20short=find(rsquare_loom_shortS20>0.3 & rsquare_loom_shortS20<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
%proportion_s20=length(idx_rsq_test_s20short)/length(rsquare_loom_shortS20);

% figure;imagesc(ZS_s20(idx_rsq_test_s20short,:), [-0.5 4]);colormap hot

end




%%

idx_rsq_test_s20short_all={};
for j=1:size(rawregressS20)
    temp_rsq=[ModelResults_shortS20_all{1,j}.rsquared];
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
    
idx_rsq_test_s20short_all{j}=idx_rsq_temp;

end

idx_rsq_test_s20short=horzcat(idx_rsq_test_s20short_all{:});
idx_rsq_test_s20short=unique(idx_rsq_test_s20short);


figure;imagesc(ZS_s20(idx_rsq_test_s20short,:), [-0.5 4]);colormap hot


%%

%%% to clasify them with a correlation

Correlation_group_s20={};
counter=1;
for i=1:size(rawregressS20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s20(idx_rsq_test_s20short),2)
        temp_corr=corrcoef(rawregressS20(i,:),ZS_s20(idx_rsq_test_s20short(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_s20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_s20=[];
for n=1:size(rawregressS20,1)
Correlation_group_mat_s20(n,:)=Correlation_group_s20{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_s20=zeros(length(Correlation_group_mat_s20),1);
for i=1:length(Correlation_group_mat_s20)
    [~,I]=max(Correlation_group_mat_s20(:,i));
    High_corr_Nb_s20(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    
    idx_temp=find(High_corr_Nb_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_rsq_test_s20short(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_rsq_test_s20short(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_rsq_test_s20short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_rsq_test_s20short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_S20_CN_each','-dpdf','-bestfit');

%%

High_corr_Nb_s20_short=High_corr_Nb_s20;
High_corr_Nb_s20_short(find(High_corr_Nb_s20==2))=1;
High_corr_Nb_s20_short(find(High_corr_Nb_s20==4))=1;

gooodmaps=[1 3 5 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_s20_short==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_rsq_test_s20short(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_rsq_test_s20short(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_rsq_test_s20short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_rsq_test_s20short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_S20_CN_each','-dpdf','-bestfit');

%%

 save('final_S20_step1.mat','-v7.3');


