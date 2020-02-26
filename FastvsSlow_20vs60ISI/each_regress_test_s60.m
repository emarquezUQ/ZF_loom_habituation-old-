
ModelResults_shortS60_all={};
for j=1:length(rawregressS20)
    
    
ModelResults_shortS60=[];
parfor i=1:size(ZS_s60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressS20(j,:)',ZS_s60(i,ZS_short_S60));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS60(i).coef=mdl.Coefficients;
    ModelResults_shortS60(i).MSE=mdl.MSE;
    ModelResults_shortS60(i).Fitted=mdl.Fitted;
    ModelResults_shortS60(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortS60_all{j}=ModelResults_shortS60;

rsquare_loom_shortS60=[ModelResults_shortS60.rsquared];


figure;histogram(rsquare_loom_shortS60,'Normalization','probability');


idx_rsq_test_s60short=find(rsquare_loom_shortS60>0.3 & rsquare_loom_shortS60<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
proportion_s60=length(idx_rsq_test_s60short)/length(rsquare_loom_shortS60);

figure;imagesc(ZS_s60(idx_rsq_test_s60short,ZS_short_S60), [-0.5 4]);colormap hot

end




%%

idx_rsq_test_s60short_all={};
for j=1:size(rawregressS20)
    temp_rsq=[ModelResults_shortS60_all{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.5 & temp_rsq<1);
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
    
idx_rsq_test_s60short_all{j}=idx_rsq_temp;

end

idx_rsq_test_s60short=horzcat(idx_rsq_test_s60short_all{:});
idx_rsq_test_s60short=unique(idx_rsq_test_s60short);


figure;imagesc(ZS_s60(idx_rsq_test_s60short,ZS_short_S60), [-0.5 4]);colormap hot


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
    
    idx_temp=find(High_corr_Nb_s60==i);
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_rsq_test_s60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_rsq_test_s60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%%

%%% if i try to clasify with the highest r2 value.... it doesnt work as
%%% good as the correlation, specially for the inhibition regressor. 


High_rsq_Nb=zeros(size(idx_rsq_test_s60short));
for i=1:size(ZS_s60(idx_rsq_test_s60short),2)
 
 High_rsq_temp=[];   
for j=1:size(rawregressS20,1)
    
temp_rsq=[ModelResults_shortS60_all{1,j}(1,idx_rsq_test_s60short(i)).rsquared];
High_rsq_temp(1,j)= temp_rsq;

end

[~,max_idx]=max(High_rsq_temp);
High_rsq_Nb(i)=max_idx;

end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    
    idx_temp=find(High_rsq_Nb==i);
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_rsq_test_s60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_rsq_test_s60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end
