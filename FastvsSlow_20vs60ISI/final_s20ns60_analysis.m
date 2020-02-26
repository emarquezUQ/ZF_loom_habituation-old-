%%% this script is to do the possible final analysis in s20 and s60 fish. I
%%% will use the s20 denoised regressors plus the inhibition regressor
%%% found in the 200CL cluster. but I will be using the noised data for the
%%% analysis. 

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab


%%% First I load the relevant data.

%load('s20_postKmeans_CN.mat');

% this is for the inhibition regressor

load('s20_Kmeans_ZS_CN_corr010_200CL.mat');


inhib_s20_regress=Cmap_ZS_CN_corr010_200CL(110,:);

figure;plot(inhib_s20_regress);

save('inhib_s20_regress_200CL.mat','inhib_s20_regress');


%%% now the ZS with the noise, matfiles,regressors, planes and fish. i think at the
%%% moment i am not excluding any of the Slow loom fish...
load('inhib_s20_regress_200CL.mat','inhib_s20_regress');


S60data_CN=load('s60_postKmeans_CN.mat','ZS_CN','MatFiles');
ZS_s60=S60data_CN.('ZS_CN');
MatFiles_s60=S60data_CN.('MatFiles');
% ZS_s60_CN=S60data_CN.('ZS_CN');

S60data=load('s60_r2050_CL3.mat','ZS','MatFiles','idx_Plane','idx_Fish');

idx_Plane_s60=S60data.('idx_Plane');
idx_Fish_s60=S60data.('idx_Fish');

% ZS_s60=S60data.('ZS');
% MatFiles_s60=S60data.('MatFiles');


S20data_CN=load('s20_postKmeans_CN.mat','ZS_CN','MatFiles');
ZS_s20=S20data_CN.('ZS_CN');
MatFiles_s20=S20data_CN.('MatFiles');

S20data=load('s20_r2050_CL6.mat','ZS','idx_Plane','idx_Fish');

idx_Plane_s20=S20data.('idx_Plane');
idx_Fish_s20=S20data.('idx_Fish');

%ZS_s20=S20data.('ZS');

rawregressS20=load('rawregressS20.mat','rawregress');
rawregressS20 = rawregressS20.('rawregress');

rawregressS60=load('rawregressS60.mat','rawregress');
rawregressS60 = rawregressS60.('rawregress');

%%

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

%figure;plot(rawregressS60(1,ZS_short_S60)); %%% to check


%%

%%% this is for a linear regression using the s20 denoised data regressors
%%% from the 50CL + the inhibition cluster

counter=1;
figure;
for i=1:length(rawregressS20)
subplot(4,2,counter);plot(rawregressS20(i,:));
counter=counter+1;
end

rawregressS20(7,:)=inhib_s20_regress;




ModelResults_shortS60=[];
parfor i=1:size(ZS_s60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressS20',ZS_s60(i,ZS_short_S60));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS60(i).coef=mdl.Coefficients;
    ModelResults_shortS60(i).MSE=mdl.MSE;
    ModelResults_shortS60(i).Fitted=mdl.Fitted;
    ModelResults_shortS60(i).rsquared=mdl.Rsquared.Adjusted;
     
end

rsquare_loom_shortS60=[ModelResults_shortS60.rsquared];


figure;histogram(rsquare_loom_shortS60,'Normalization','probability');


idx_rsq_test_s60short=find(rsquare_loom_shortS60>0.3 & rsquare_loom_shortS60<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
proportion_s60=length(idx_rsq_test_s60short)/length(rsquare_loom_shortS60);

figure;imagesc(ZS_s60(idx_rsq_test_s60short,ZS_short_S60), [-0.5 4]);colormap hot



%%


coefficients_s60={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(ModelResults_shortS60)%%% to make a variable the size of ModelResults
    coef=[ModelResults_shortS60(idx).coef];%%%% and then put in another variable the coef field from ModelResults
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    if ~isempty(temp)%%% if temp is not empty...
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            if coef.pValue(coef_idx)<0.05%%%to select the coef that are bellow the p value we want, in this case 0.05
                coefficients_s60{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients_s60); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
coefficients_s60(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
coefficients_s60=cell2mat(coefficients_s60); %%%to make a matrix of the coefficients array


coefficients_Nb_s60=coefficients_s60>0; %%%to take the values above 0 and make a variable
coefficients_Nb_s60=sum(coefficients_Nb_s60,2);%%% to sum the cells from the second dimension (the looms)
idx_coef_s60=find(coefficients_Nb_s60>0);%%%to take that have more than 2 responses (that how it was before) but i might want to try to see what happens if i dont do this
idx_coef_rsq_s60=intersect(idx_rsq_test_s60short,idx_coef_s60);%%% to make a variable with the common values from idx_rsq and idx_coef, so with the r2 values between 0.3-1 and the coefficients that were significant



%%

%%% to clasify them with a correlation

Correlation_group_s60={};
counter=1;
for i=1:size(rawregressS20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s60(idx_coef_rsq_s60),1)
        temp_corr=corrcoef(rawregressS20(i,:),ZS_s60(idx_coef_rsq_s60(idx),ZS_short_S60));
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
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_coef_rsq_s60(idx_temp),ZS_short_S60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_coef_rsq_s60(idx_temp),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_coef_rsq_s60(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_coef_rsq_s60(idx_temp)));%%% for the fish location
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
    
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_newS20(i).coef=mdl.Coefficients;
    ModelResults_newS20(i).MSE=mdl.MSE;
    ModelResults_newS20(i).Fitted=mdl.Fitted;
    ModelResults_newS20(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom_newS20=[ModelResults_newS20.rsquared];
figure;histogram(rsquare_loom_newS20,'Normalization','probability');

idx_rsq_test_new_s20=find(rsquare_loom_newS20>0.3 & rsquare_loom_newS20<1);%%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
proportion_s20=length(idx_rsq_test_new_s20)/length(rsquare_loom_newS20)

figure;imagesc(ZS_s20(idx_rsq_test_new_s20,:), [-0.5 4]);colormap hot

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
    
    idx_temp=find(High_corr_Nb_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_rsq_test_new_s20(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_rsq_test_new_s20(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_rsq_test_new_s20(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_rsq_test_new_s20(idx_temp)));%%% for the fish location
    counter=counter+4;
end



%  save('final_S20nS60_step1.mat','-v7.3');


%%

%%% now I am going to try to get the idx of the ROIs depending of how much
%%% they fit the regressors and build groups with that. 


%%% I first tried with the hights coef. before didnt work well...
%%% then i will triy with the most significant as before it didnt work very well
%%% either, although a bit better.
%%% 
High_coeff_s20={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(ModelResults_newS20(idx_rsq_test_new_s20))%%% to make a variable the size of ModelResults
    coef=[ModelResults_newS20(idx_rsq_test_new_s20(idx)).coef];%%%% and then put in another variable the coef field from ModelResults
    %coef2=[ModelResults(idx_coef_rsq(idx)).coef.Estimate];
    %coef2=sortrows(coef2(2:height(coef)),'descend');
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    
    %min_pv=min(table2array(coef(2:height(coef),4)));
    max_est=max(table2array(coef(:,1)));
    if ~isempty(temp)%%% if temp is not empty...
        
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            if coef.Estimate(coef_idx)==max_est && coef.Estimate(coef_idx)
                %coef.pValue(coef_idx)==min_pv
                %coef.Estimate(coef_idx)==max_est
                High_coeff_s20{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            else  
            
            end
        end
    end
end
idxempty=cellfun('isempty',High_coeff_s20); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
High_coeff_s20(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
High_coeff_s20=cell2mat(High_coeff_s20); %%%to make a matrix of the coefficients array


%High_coeff_Nb=High_coeff>0; %%%to take the values above 0 and make a variable
High_coeff_Nb_s20=zeros(length(High_coeff_s20),1);
for i=1:length(High_coeff_s20)
    High_coeff_Nb_s20(i,1)=find(High_coeff_s20(i,:),1,'first');
    
   %High_coeff_Nb(i,1)=subsref(find(High_coeff(i,:)),struct('type','()','subs',{{1}}));
end


idxGroup_final_s20=zeros(length(ModelResults_newS20),1);
idxGroup_final_s20(idx_rsq_test_new_s20)=High_coeff_Nb_s20;



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    
    idx_temp=find(High_coeff_Nb_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_rsq_test_new_s20(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_rsq_test_new_s20(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_rsq_test_new_s20(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_rsq_test_new_s20(idx_temp)));%%% for the fish location
    counter=counter+4;
end

