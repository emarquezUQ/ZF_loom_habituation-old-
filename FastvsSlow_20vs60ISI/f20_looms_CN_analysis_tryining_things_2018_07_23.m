%%%% this is to try a PCA in the FnS_20vs60ISI data, using also the Noise...
%%%% first triying with the f20

MatFiles=dir('*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium
%MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:);  %%%to combine the Calcium and Fitness variables (need to ask Gilles what Fitness is). Fitness here is the variable were we take the good calcium responses from the HPC analysis and pairthem with their index number.
%GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);


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

% ZS=zscore(GoodCalcium,1,2); %%%to normalize the data
% ZS=detrend(ZS')';%%% to Remove a linear trend from ZS (why???)

%%% this was to make a PCA, but we have checked the results and is not that
%%% interesting
% [coeff_ZS,score_ZS,~,~,explained_ZS,~] = pca(ZS);
% figure;plot(cumsum(explained_ZS)); 

ZS_CN=zscore(GoodCalNoise,1,2); %%%to normalize the data
ZS_CN=detrend(ZS_CN')';%%% to Remove a linear trend from ZS (why???)

%%% this was to make a PCA, but we have checked the results and is not that
%%% interesting
% [coeff,score,~,~,explained,~] = pca(ZS_CN);
% figure;plot(cumsum(explained)); 


 clear Calcium Fitness GoodCalcium Noise GoodCalNoise
% 
%  save('f20_CN.mat','-v7.3');
 
%%
%%%here i am doing a kmeans to get some raw regressors
options = statset('UseParallel',1); [idxKmeans_ZS_CN Cmap_ZS_CN]=kmeans(ZS_CN,50,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS_CN2,GoodBetas_ZS_CN2]=Test_Regress(Cmap_ZS_CN,Cmap_ZS_CN,idxKmeans_ZS_CN,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

saveas(gcf,'50rawclusters_CN_f20','jpg');


% clear Calcium Fitness GoodCalcium
% 
%  save('f20_postKmeans_CN.mat','-v7.3');
 %%
 goodmaps_CN=[7 22 29 37 45 47];
for i=1:length(goodmaps_CN)
rawregress_CN(i,:)=Cmap_ZS_CN(goodmaps_CN(i),:);
end

%%% i am adding a multisensory regressor by adding part of a fasthab
%%% cluster and the soundspondent cluster (37 and 47). 

multisensregess=zeros(1,size(rawregress_CN,2));
multisensregess(1,1:140)=rawregress_CN(5,1:140)/2;
multisensregess(1,141:1344)=rawregress_CN(4,141:1344);

figure;plot(multisensregess);
rawregress_CN(7,:)=multisensregess;

multisensregess2=zeros(1,size(rawregress_CN,2));
multisensregess2(1,:)=rawregress_CN(3,:)*2;
multisensregess2(1,900:1000)=rawregress_CN(3,900:1000)*2+rawregress_CN(4,900:1000);

figure;plot(rawregress_CN(4,:));
figure;plot(multisensregess2);
rawregress_CN(8,:)=multisensregess2;

save ('rawregressF20_CN.mat','rawregress_CN');

% for i=1:length(goodmaps)
% figure;
% plot(rawregress(i,:));
% end


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(goodmaps_CN);

for i=1:length(goodmaps_CN)+2 %%%for the multisensregress
   
    subplot(4,2,counter);plot(rawregress_CN(i,:));
     counter=counter+1;
end

saveas(gcf,'goodmaps_F20_CN_3','jpg');


%%%now i will try to look if some of the regressors are only in one fish to
%%%select the best ones. 

Numbers=[1 [MatFiles.GoodNumber]]; %%%to take make a vector with the GoodNumber field
counter=1;
idx_Plane=nan(length(ZS_CN),1);%%%% to make an empty (with nans) one column variable the size of ZS
idx_Fish=nan(length(ZS_CN),1);%%%% to make an empty (with nans) one column variable the size of ZS
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %%%to get the number of the plane    
    idx_Plane(Numbers(i):Numbers(i+1))=Plane; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Plane   
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1}); %%%to get the number of the fish 
    idx_Fish(Numbers(i):Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore


%%%this is to make a figure where we will plot the mean of clusters
%%%selected the raster plot, the levels and the fish where they are found.
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(goodmaps_CN);counter=1;
for i=goodmaps_CN
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans_ZS_CN==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_temp,:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_temp));%%% for the fish location
    counter=counter+4;
end
saveas(gcf,'raw_multigraph_F20_CN','jpg');


ModelResults=[];
parfor i=1:size(ZS_CN,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregress_CN',ZS_CN(i,:));
    %mdl=stepwiselm(rawregress',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    %%%this is to put the results in the ModelResulsts variable
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them
idx_rsq=find(rsquare_loom>0.5 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; %%%to plot them in a raster plot
imagesc(ZS_CN(idx_rsq,:), [-0.5 4]);colormap hot

saveas(gcf,'rasterplot_r050_f20_CN_3','jpg');

figure; 
plot(mean(ZS_CN(idx_rsq,:)));




ZS_CNrsq05=ZS_CN(idx_rsq,:);

save('f20_postRegress_CN_R2_050.mat','-v7.3');

%%


coefficients={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(ModelResults)%%% to make a variable the size of ModelResults
    coef=[ModelResults(idx).coef];%%%% and then put in another variable the coef field from ModelResults
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    if ~isempty(temp)%%% if temp is not empty...
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            if coef.pValue(coef_idx)<0.05%%%to select the coef that are bellow the p value we want, in this case 0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
coefficients(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
coefficients=cell2mat(coefficients); %%%to make a matrix of the coefficients array


coefficients_Nb=coefficients>0; %%%to take the values above 0 and make a variable
coefficients_Nb=sum(coefficients_Nb,2);%%% to sum the cells from the second dimension (the looms)
idx_coef=find(coefficients_Nb>0);%%%to take that have more than 2 responses (that how it was before) but i might want to try to see what happens if i dont do this
idx_coef_rsq=intersect(idx_rsq,idx_coef);%%% to make a variable with the common values from idx_rsq and idx_coef, so with the r2 values between 0.3-1 and the coefficients that were significant

%%
%%% this was to make a PCA, but we have checked the results and is not that
%%% interesting
% [coeff_rsq,score_rsq,~,~,explained_rsq,~] = pca(ZS_CN(idx_coef_rsq,:));
% figure;plot(cumsum(explained_rsq)); 



% options = statset('UseParallel',1); [idxKmeans_ZS_CN_rsq Cmap_ZS_CN_rsq]=kmeans(ZS_CN(idx_coef_rsq,:),20,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
% [Model_ZS_CN2,GoodBetas_ZS_CN2]=Test_Regress(Cmap_ZS_CN_rsq,rawregress_CN,idxKmeans_ZS_CN_rsq,0.5);%%%here we do another linear regression and we select the ones with an r2 value above 0.3
%  
% saveas(gcf,'goodclusters_CN_r050_f20_CL7','jpg');

%%% when including the noise data the Kmeans after the regression takes
%%% away the "inhibition" cluster... so in this case I wont use the 2nd
%%% kmeans clusters and and will keep the original ones. 

%%% Note:but now I need to add also the multisens ROIs


 

idxKmeans1_coef_rsq_CN=idxKmeans_ZS_CN(idx_coef_rsq);%%% to make a variable with the common idx of the interesting clusters from the 1st kmeans and the traces hat passed the regression
rows=length(goodmaps_CN);
counter=1;

%%%this is to make a figure where we will plot the mean of clusters
%%%selected the raster plot, the levels and the fish where they are found.
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(goodmaps_CN);counter=1;
for i=goodmaps_CN
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq_CN==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_coef_rsq(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end
saveas(gcf,'multigraph_CN_r050_F20','jpg');

Test=[];i=1;
name=strcat(MatFiles(i).name);
Rs=load(name, 'ROIs');%%% to load the ROIs
Rs=Rs.ROIs;%%% and put them on a variable
F=load(name, 'idx_components');%%% to load the ROIs idx_components again
F=F.idx_components+1;%%% and put them on a variable
Rs=Rs(:,F);%%% to combine the ROis and the idx_components
MatFiles(i).ROI=Rs;%%% to make an ROI field in MatFiles
Test(i)=size(Rs,2)==MatFiles(i).GoodNumber; %%% to check of the size of the columns of the ROIs is the same as the first (cause i=1) goodnumber 
for i = 2:length(MatFiles)%%% to take the planes (except the first one as we already did it... ) and do the same as above
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    MatFiles(i).ROI=Rs;
    Test(i)=size(Rs,2)==(MatFiles(i).GoodNumber-MatFiles(i-1).GoodNumber); %%%to to check of the size of the columns of the ROIs is the same as the goodnumber in turn. when they are you get a 1 in a column in Test until you get 25 columns. but i dont understand why we do it...
end
clearvars GC C S F N name i;

idxKmeans_final=zeros(size(ZS_CN,1),1);%%%to make a variable with 0s the size of ZS_CN
idxKmeans_final(idx_coef_rsq)=idxKmeans1_coef_rsq_CN;%%% and put the filtered responses by r2 and filters

temp=[];
counter=1;
for i=goodmaps_CN %%%to take the clusters we want
    idx_temp=find(idxKmeans_final==i);   %%%and put the filtered values in a new variable in different columns
    temp{counter}=idx_temp;    
    counter=counter+1;    
end

Numbers(1)=0; %%% to change the first value of Numbers to 0 (cause it was not relevant)
%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(goodmaps_CN),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
% colors = [0         0    1.0000
%          0    0.5000    1.0000
%     1.0000         0         0
%     1.0000    0.1034    0.7241
%     1.0000    0.5000    0.3000
%          0    0.7000    0.2000
%     0.5000    0.5000         0
%          0    0.5000    0.5000];
colors = colors*256;


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(goodmaps_CN)));yplot=ceil(length(GoodBetas_ZS_CN2)/xplot);
for i=goodmaps_CN
    idx_temp=find(idxKmeans1_coef_rsq_CN==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_coef_rsq(idx_temp),:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end

saveas(gcf,'good_multigraph_CN_r050_F20_2','jpg');

% clear Calcium Fitness GoodCalcium
% 
%  save('f20_CN_r2050_CL5.mat','-v7.3');
 
 %%
 %%% Gilles and I explore the PCA results but they were not really
 %%% interesting so I will give up on them for the moment. 
%  figure;
% subplot(1,3,1);imagesc(coeff_ZS(:,1:20));
% subplot(1,3,2);imagesc(coeff(:,1:20));
% subplot(1,3,3);imagesc(coeff_rsq(:,1:20));
% 
% figure;
% plot(cumsum(explained_rsq));
% 
% figure;
% max_i=20;
% for i=1:max_i
%     subplot(max_i,1,i);plot(coeff_rsq(:,i));
% end
% 
% figure;
% score_idx=2;
% plot(mean(ZSrsq05(find(score_rsq(:,score_idx)>(mean(score_rsq(:,score_idx))+1*std(score_rsq(:,score_idx)))),:),1));title(num2str(length(find(score_rsq(:,score_idx)>(mean(score_rsq(:,score_idx))+1*std(score_rsq(:,score_idx)))))));
% 
% figure;
% for score_idx=1:10;
% subplot(10,1,score_idx);
% plot(mean(ZSrsq05(find(score_rsq(:,score_idx)>(mean(score_rsq(:,score_idx))+1*std(score_rsq(:,score_idx)))),:),1));title(num2str(length(find(score_rsq(:,score_idx)>(mean(score_rsq(:,score_idx))+1*std(score_rsq(:,score_idx)))))));
% end
% 
% [L1,T] = rotatefactors(coeff_rsq(:,1:20),'Method','promax');
% figure;imagesc(L1);
% 
% 
% 
% [~,test]=Test_Regress(coeff',rawregress_CN,ones(size(coeff)),0.1);
% 
% 
% figure;plot(mean(ZSrsq05(find(idxKmeans1_coef_rsq_CN==45),:)),1);hold on;plot(mean(ZSrsq05(find(idxKmeans1_coef_rsq_CN==47),:)),1);ylim([0 5]);

 
 
 %%
%%%%%%%%%%%%%%%%%%%%%%%%%% this part is to make the Kmeans figures with the localization of the selected clusters in the brain  %%%%%%%%%
idx=1
    filename=MatFiles(idx).name;
%destdirectory = strcat('/QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab','/Kmeans_f20_r2050_CL50_HPC',filename(8:11),'/');
destdirectory = strcat('C:/Emm_temp/FvsS_20vs60_analysis','/Kmeans_f20_CN_r2050_CL50',filename(8:11),'/');
mkdir(destdirectory);   %create the directory
clearvars idx filename

for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;%%%to put the respective name of the files (in this case the slices)
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
    %for k = 1 : length(temp)
    %%%%%this is to locate in every plane which ROIs and cluster can be
    %%%%%found
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1) & [temp{k}]>Numbers(idx)); %%%to put the data filtered of the selected clusters on a new variable. but i dont understand the sentece... why not just ask for Numbers(idx)?
        if tempROIsNb            
            ROIsNb=[ROIsNb; temp{k}(tempROIsNb)];%%%%
            %temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb; repmat(k,length(tempROIsNb),1)];
        end
    end
    
    %%%% and this part is to make and image using the _mean.tif image as
    %%%% base, where we will add the ROIs located before and colored base
    %%%% on the number of the cluster
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles(idx).ROI;       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        %%%finally we save the Image
        name=strcat(destdirectory,'_Kmeans_f20_CN_r2050_CL50',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end


clearvars idx i tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

%%

%%% now I am going to try to get the idx of the ROIs depending of how much
%%% they fit the regressors and build groups with that. 


% rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them
% idx_rsq=find(rsquare_loom>0.5 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1

%%% I first tried with the hights coef but didnt work well...
%%% then i will tried with the most significant but it didnt work very well
%%% either, although a bit better.
%%% 
High_coeff={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(ModelResults(idx_coef_rsq))%%% to make a variable the size of ModelResults
    coef=[ModelResults(idx_coef_rsq(idx)).coef];%%%% and then put in another variable the coef field from ModelResults
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
                High_coeff{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            else  
            
            end
        end
    end
end
idxempty=cellfun('isempty',High_coeff); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
High_coeff(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
High_coeff=cell2mat(High_coeff); %%%to make a matrix of the coefficients array


%High_coeff_Nb=High_coeff>0; %%%to take the values above 0 and make a variable
High_coeff_Nb=zeros(length(High_coeff),1);
for i=1:length(High_coeff)
    High_coeff_Nb(i,1)=find(High_coeff(i,:),1,'first');
    
   %High_coeff_Nb(i,1)=subsref(find(High_coeff(i,:)),struct('type','()','subs',{{1}}));
end


idxGroup_final=zeros(length(ModelResults),1);
idxGroup_final(idx_coef_rsq)=High_coeff_Nb;



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregress_CN,1);counter=1;
for i=1:size(rawregress_CN,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(High_coeff_Nb==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_coef_rsq(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end


%%
%%% here I am trying to use the coefficients to sort the ROIs using Kmeans
%%% it didnt work very well...
options = statset('UseParallel',1); [idxKmeans_coeff_CN_rsq Cmap_coeff_CN_rsq]=kmeans(coefficients(idx_coef_rsq,:),7,'Options',options,'Replicates',5,'MaxIter',1000,'Display','final');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=max(idxKmeans_coeff_CN_rsq);counter=1;
for i=1:max(idxKmeans_coeff_CN_rsq)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans_coeff_CN_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_coef_rsq(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%%
%%% trying to use the classification learner app
%%% first I need to transfomr into a table my dataset

%ZS_CNrsq05_table = array2table(ZS_CNrsq05);

%%% it seems that the best method for our set of data is a K-NN. 

%%% I could use 2 sets of data, the coefficient matrix OR the ZS_rsq matrix
%%% -the coefficient matrix would be easier and quicker but might not give
%%% the best results. This one would work better with euclidean distance.
%%%-the ZS_rsq matrix will give better results but it could crash teh Iron
%%%Horse and HPC jobs... This one would work better with cityblock distance.

%%% either way I need a matrix or table with a sample (e.g.20-30ROIs) and
%%% its values (ZS or coefficients) AND a vector with the classification of
%%% those ROIs.
%%%(I chose that classification somehow, i will try after a very strong correlation or
%%% take the idx from the kmeans.
%%
%%% Here I am trying to look for the specific neuron groups using a linear
%%% regressin... Is not working well neither with the inhibition group or the multisens group...


ModelResults_temp=[];
parfor i=1:size(ZS_CN,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregress_CN(7,:)',ZS_CN(i,:));
    %mdl=stepwiselm(rawregress',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_temp(i).coef=mdl.Coefficients;
    ModelResults_temp(i).MSE=mdl.MSE;
    ModelResults_temp(i).Fitted=mdl.Fitted;
    ModelResults_temp(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom=[ModelResults_temp.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them
idx_rsq=find(rsquare_loom>0.62 & rsquare_loom<1); %%%then select the rsquare that are between 0.3 and 1
figure; %%%to plot them in a raster plot
imagesc(ZS_CN(idx_rsq,:), [-0.5 4]);colormap hot

figure; %%%to plot them in a raster plot
plot(mean(ZS_CN(idx_rsq,:)));

%%
%%% and this is to try the same using just a correlation. 
%%% is working OKish for the inhibition group but not great for multisens.
%%% after modifiying a bit the multisens regressor is working better. 0.75
%%% gives good results but very small number of ROIs. 0.5 gives some Okish
%%% resulsts too. 

 
Correlation_group={};
counter=1;
for i=1:size(rawregress_CN,1)
    Correlation=[];
    for idx=1:size(ZS_CNrsq05,1)
        temp_corr=corrcoef(rawregress_CN(i,:),ZS_CNrsq05(idx,:));
            Correlation(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group{i}=Correlation;
    counter=counter+1;
end

Correlation_group_mat=[];
for n=1:size(rawregress_CN,1)
Correlation_group_mat(n,:)=Correlation_group{n}(n,:);

end

Correlation_group_all={};
counter=1;
for i=1:size(rawregress_CN,1)
    Correlation=[];
    for idx=1:size(ZS_CN,1)
        temp_corr=corrcoef(rawregress_CN(i,:),ZS_CN(idx,:));
            Correlation(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_all{i}=Correlation;
    counter=counter+1;
end

clear Correlation

for n=1:length(Correlation_group_all)
Correlation_group_all_mat(n,:)=Correlation_group_all{n}(n,:);

end



%%%%% this is to plot 

Threshold=0.5;
   
    for n=1:length(Correlation_group_all)
         figure;
    
        Correlation=Correlation_group{n};
        plot(mean(ZS_CNrsq05(find(Correlation(n,:)>Threshold),:)));
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        mdl=fitlm(rawregress_CN(n,:)',mean(ZS_CNrsq05(find(Correlation(n,:)>Threshold),:))); %%% I changed it to fitlm cause is faster
        title(strcat('lin_reg to looms : ',num2str([mdl.Rsquared.Adjusted]),' nb of ROIs : ', num2str(length(find(Correlation(n,:)>Threshold)))));
       figure; %%%to plot them in a raster plot
       imagesc(ZS_CNrsq05(find(Correlation(n,:)>Threshold),:), [-0.5 4]);colormap hot
    end
    

    

    
%%
%%% now I will try to do a Kmeans on the correlation matrix and see if it
%%% works

%%% it doesnt work that bad... the best so far, there was just a problem with the multisens
%%% group but the rest is actually not that bad


options = statset('UseParallel',1); [idxKmeans_corr_CN_rsq Cmap_corr_CN_rsq]=kmeans(Correlation_group_mat',7,'Options',options,'Replicates',5,'MaxIter',1000,'Display','final');


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=max(idxKmeans_corr_CN_rsq);counter=1;
for i=1:max(idxKmeans_corr_CN_rsq)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans_corr_CN_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_rsq(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end

saveas(gcf,'multigraph_CN_corr_F20','jpg');



%%%what if I filter looking for the max correlation?
%%% it seems to be working!!!!

High_corr_Nb=zeros(length(Correlation_group_mat),1);
for i=1:length(Correlation_group_mat)
    [~,I]=max(Correlation_group_mat(:,i));
    High_corr_Nb(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregress_CN,1);counter=1;
for i=1:size(rawregress_CN,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(High_corr_Nb==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_rsq(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end

saveas(gcf,'multigraph_CN_maxcorr_F20_3','jpg');

  save('f20_CN_r2050_CL5_extra.mat','-v7.3');
  
  %%
%%% now I will try the correlation in the whole data set to and filter them
%%% by 0.5 to see if the number of ROIs change compare to filtering by the
%%% linear regression (see above).
%%% it works fine it seems, is more permisive but almost increases the
%%% number of ROIs 3 times. 

 
Correlation_group_all={};
counter=1;
for i=1:size(rawregress_CN,1)
    Correlation=[];
    for idx=1:size(ZS_CN,1)
        temp_corr=corrcoef(rawregress_CN(i,:),ZS_CN(idx,:));
            Correlation(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_all{i}=Correlation;
    counter=counter+1;
end

clear Correlation

for n=1:size(rawregress_CN,1)
Correlation_group_all_mat(n,:)=Correlation_group_all{n}(n,:);

end

%%% first I need to pick the ones higher than 0.5

temp_idx_corr=zeros(length(Correlation_group_all_mat),1);
for i=1:length(Correlation_group_all_mat)
    if length(find(Correlation_group_all_mat(:,i)>0.5))>0
     
    temp_idx_corr(i,1)=1;
    else
    temp_idx_corr(i,1)=0;    
    end    
end

idx_corr=find(temp_idx_corr);

%%% testing with higher than 0.75... is too much, I get less ROIs that the
%%% linear regression (with 0.5).  a correlation with 0.7 gives a similar result than the 0.5
%%% linear regression.

temp_idx_corr=zeros(length(Correlation_group_all_mat),1);
for i=1:length(Correlation_group_all_mat)
    if length(find(Correlation_group_all_mat(:,i)>0.7))>0
     
    temp_idx_corr(i,1)=1;
    else
    temp_idx_corr(i,1)=0;    
    end    
end

idx_corr070=find(temp_idx_corr);

%%%now I filter looking for the max correlation


High_corr_Nb_all=zeros(length(Correlation_group_all_mat(idx_corr)),1);
for i=1:length(Correlation_group_all_mat(idx_corr))
    [~,I]=max(Correlation_group_all_mat(:,idx_corr(i,1)));
    High_corr_Nb_all(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregress_CN,1);counter=1;
for i=1:size(rawregress_CN,1)
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(High_corr_Nb_all==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_corr(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_corr(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_corr(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_corr(idx_temp)));%%% for the fish location
    counter=counter+4;
end

saveas(gcf,'multigraph_CN_maxcorr_F20_all_3','jpg');


%%
%%% making new regressors to find incrase in activity (possible gabaergic
%%% neurons inhibiting other ones). 


%%% for f20. first we set the times
increaseregress1=zeros(6720,1);
startpoints=[30,254,478]; %%in seconds
loom_times=[0,22,40,60,78,100,120,140,158,180]; %%in seconds
for p=1:3
for k=1:10
    increaseregress1(startpoints(p)*10+loom_times(k)*10)=2;
end
end
temp_incr1=find(increaseregress1);
temp_incr1=round(temp_incr1/5);

%the we put the spikes.
increaseregress1=zeros(1,1344);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';

for i=1:30
    increaseregress1(temp_incr1(i):temp_incr1(i)+size(GCaMP6,1)-1)=GCaMP6/6;
end


quickincreaseregress=increaseregress1;
%%% and i correct first for the quick increase
p=[96 24 12 12 12];
counter=1;
for i=[1 2 3 11 21]
    quickincreaseregress(temp_incr1(i):temp_incr1(i)+size(GCaMP6,1)-1)=GCaMP6/(p(counter));
    counter=counter+1;
end

figure; plot(quickincreaseregress);


slopeincreaseregress=increaseregress1;
%%% and i correct first for the quick increase
counter=1;
for i=4
    slopeincreaseregress(temp_incr1(i):temp_incr1(i)+size(GCaMP6,1)-1)=(GCaMP6/(6-counter));
    slopeincreaseregress(temp_incr1(i+10):temp_incr1(i+10)+size(GCaMP6,1)-1)=(GCaMP6/(6-0.5));
    slopeincreaseregress(temp_incr1(i+20):temp_incr1(i+20)+size(GCaMP6,1)-1)=(GCaMP6/(6-0.5));
    counter=counter+0.2;
end

figure; plot(slopeincreaseregress);

increaseregress(1,:)=quickincreaseregress;
increaseregress(2,:)=slopeincreaseregress;

%%% now lets try to find ROIs that correlate to them

%%% I dont find many neurons with increase of activity (less that 20 if Threshold is 0.7 and less than 200 with threshold of 0.6).
%%% if i reduce the threshcold i get just non-hab neurons...

%%%note... i moved something and now is not working properlly... i think i
%%%found the problem. is how i am indexing in the Correlation part... is
%%%not the same as before... compare Correlation_group_in1 vs 2.
%%%solved the problem! it was becuase for the Correlation_group_in1 i was
%%%using all the ROIs and for the Correlation_group_in2 i was using only
%%%the ones that correlate to non-hab
 

Correlation_group_in1={};
counter=1;

for i=1:2
    Correlation=[];
    for x=1:size(ZS_CN,1)%%%length(idx_nonhab_corr)
        %idx=idx_nonhab_corr(x);
        temp_corr=corrcoef(increaseregress(i,:),ZS_CN(x,:));
            Correlation(i,x)=temp_corr(1,2);
                        
        
    end
    Correlation_group_in1{counter}=Correlation;
    counter=counter+1;
end



idx_nonhab_corr=idx_corr(find(High_corr_Nb_all==3));

figure; plot((mean(ZS_CN(idx_nonhab_corr,:))));

Correlation_group_in2={};
counter=1;

for i=1:2
    Correlation=[];
    for x=1:length(idx_nonhab_corr)
        %idx=idx_nonhab_corr(x);
        temp_corr=corrcoef(increaseregress(i,:),ZS_CN(idx_nonhab_corr(x),:));
            Correlation(i,x)=temp_corr(1,2);
                        
        
    end
    Correlation_group_in2{counter}=Correlation;
    counter=counter+1;
end


%figure;plot(mean(ZS_CN(idx_corr(find(High_corr_Nb_all==3)),:),1));
%figure;plot(mean(ZS_CN(idx_corr(find(High_corr_Nb_all==3)),1:200),1));


%%%%% this is to plot 

Threshold=0.5;
   
    for n=1:2
         figure;
    
        Correlation=Correlation_group_in2{n};
        plot(mean(ZS_CN(find(Correlation(n,:)>Threshold),:)));
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        mdl=fitlm(increaseregress(n,:)',mean(ZS_CN(find(Correlation(n,:)>Threshold),:))); %%% I changed it to fitlm cause is faster
        title(strcat('lin_reg to looms : ',num2str([mdl.Rsquared.Adjusted]),' nb of ROIs : ', num2str(length(find(Correlation(n,:)>Threshold)))));
       figure; %%%to plot them in a raster plot
       imagesc(ZS_CN(find(Correlation(n,:)>Threshold),:), [-0.5 4]);colormap hot
    end
    
    
    Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(Correlation_group_in2);counter=1;
for i=1:size(Correlation_group_in2,2)
    Correlation=Correlation_group_in2{i};
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(Correlation(i,:)>Threshold);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_temp,:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_temp));%%% for the fish location
    counter=counter+4;
end

 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%% this part is to make the Kmeans figures with the localization of the selected clusters in the brain  %%%%%%%%%

temp=[];
counter=1;
for i=1:size(Correlation_group_in2,2)
    Correlation=Correlation_group_in2{i};
temp{counter}=find(Correlation(i,:)>Threshold);    
counter=counter+1;    

end

Numbers(1)=0; %%% to change the first value of Numbers to 0 (cause it was not relevant)
%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(temp),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
% colors = [0         0    1.0000
%          0    0.5000    1.0000
%     1.0000         0         0
%     1.0000    0.1034    0.7241
%     1.0000    0.5000    0.3000
%          0    0.7000    0.2000
%     0.5000    0.5000         0
%          0    0.5000    0.5000];
colors = colors*256;


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(temp)));yplot=ceil(length(temp)/xplot); rows=length(temp);
for i=1:length(temp)
    idx_temp=temp{i};
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_temp,:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_temp));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end

%saveas(gcf,'multigraph_LoomMov_Mov_CN_F20_R04','jpg');

idx=1
    filename=MatFiles(idx).name;
%destdirectory = strcat('/QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab','/Kmeans_f20_r2050_CL50_HPC',filename(8:11),'/');
destdirectory = strcat('C:/Emm_temp/FvsS_20vs60_analysis','/Kmeans_f20_CN_increase_R06',filename(8:11),'/');
mkdir(destdirectory);   %create the directory
clearvars idx filename

for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;%%%to put the respective name of the files (in this case the slices)
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
    %for k = 1 : length(temp)
    %%%%%this is to locate in every plane which ROIs and cluster can be
    %%%%%found
    for k = length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1) & [temp{k}]>Numbers(idx)); %%%to put the data filtered of the selected clusters on a new variable. but i dont understand the sentece... why not just ask for Numbers(idx)?
        if tempROIsNb            
            ROIsNb=[ROIsNb; temp{k}(tempROIsNb)];%%%%
            %temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb; repmat(k,length(tempROIsNb),1)];
        end
    end
    
    %%%% and this part is to make and image using the _mean.tif image as
    %%%% base, where we will add the ROIs located before and colored base
    %%%% on the number of the cluster
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles(idx).ROI;       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        %%%finally we save the Image
        name=strcat(destdirectory,'_Kmeans_f20_CN_increase_R06',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end

%%
%%% here i am trying to look again form the increasing activity neurons but
%%% by find the ones that didnt respond to the first loom. 
%%% it doesnt seem to be very effective...

 idx_nonhab_corr=idx_corr(find(High_corr_Nb_all==3));

figure; plot((mean(ZS_CN(idx_nonhab_corr,:))));

size(idx_nonhab_corr,1)

Correlation_group_in2={};
counter=1;

for i=1:2
    Correlation=[];
    for x=1:size(ZS_CN,1)%%%length(idx_nonhab_corr)
        %idx=idx_nonhab_corr(x);
        temp_corr=corrcoef(increaseregress(i,:),ZS_CN(x,:));
            Correlation(i,x)=temp_corr(1,2);
                        
        
    end
    Correlation_group_in2{counter}=Correlation;
    counter=counter+1;
end


%figure;plot(mean(ZS_CN(idx_corr(find(High_corr_Nb_all==3)),:),1));
%figure;plot(mean(ZS_CN(idx_corr(find(High_corr_Nb_all==3)),1:200),1));

figure;
plot(mean(ZS_CN(find(Correlation(n,:)>Threshold),:)));



%%%%% this is to plot 

Threshold=0.6;
   
 Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(Correlation_group_in);counter=1;
for i=1:size(Correlation_group_in,2)
    Correlation=Correlation_group_in{i};
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(Correlation(i,:)>Threshold);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_temp,:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_temp));%%% for the fish location
    counter=counter+4;
end  
  
Threshold=0.3;
Correlation=Correlation_group_in{2};
idx_corr_group_in=find(Correlation(2,:)>Threshold);

figure;plot(mean(ZS_CN(idx_corr_group_in,:))); 

increase3=zeros(length(idx_corr_group_in),1);
%idx=1:length(idx_corr(find(High_corr_Nb_all==3)));
for i=1:length(idx_corr_group_in)
    if true(ZS_CN(idx_corr_group_in(i),1:90)<1)
    
    increase3(i,1)=1;
    else
        increase3(i,1)=0;
    end
end
  idx_corr_increase3=find(increase3);

   figure;plot(mean(ZS_CN(idx_corr_group_in(idx_corr_increase3),:)));  
   
    Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=1; counter=1;
%for i=1:size(Correlation_group_in,2)
    
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=idx_corr_group_in(idx_corr_increase3);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_temp,:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_temp));%%% for the fish location
    counter=counter+4;
%end  

 %%
 %%%here i am doing a kmeans to get some raw regressors but i tried to
%%%overcluster so see if in that way I get the multisensitive ones or de
%%%increasing ones. it didnt seem to work very well...

options = statset('UseParallel',1);
[idxKmeans_ZS_CN_500CL Cmap_ZS_CN_500CL]=kmeans(ZS_CN,500,'Options',options,'Distance','cityblock','Replicates',2,'MaxIter',1000,'Display','final');


%[idxKmeans_ZS_CN_500CL Cmap_ZS_CN_500CL]=kmeans(ZS_CN,500,'Distance','cityblock','Replicates',2,'MaxIter',1000,'Display','final');
%[Model_ZS_CN_500CL,GoodBetas_ZS_CN2_500CL]=Test_Regress(Cmap_ZS_CN,Cmap_ZS_CN,idxKmeans_ZS_CN_500CL,0.1);%%%here we do another linear regression and we select the ones with an r2 value above 0.3


%saveas(gcf,'500rawclusters_CN_f20','jpg');

counter2=0;
for i=1:10
    Fighandle=figure;counter=1;
    set(Fighandle, 'Position', [100, 100, 1500, 1000]);
    for j=1:50
        subplot(7,8,counter);plot(Cmap_ZS_CN_500CL(j+counter2,:));
        counter=counter+1;
    end
    counter2=counter2+50;
    print(Fighandle,strcat('500rawclusters_CN_f20_',num2str(i)),'-dpng','-r0');
    %close all;
end


save ('f20_idxKmeans_n_Cmap_ZS_CN_500CL.mat','idxKmeans_ZS_CN_500CL','Cmap_ZS_CN_500CL');

%%

%%% Here I filtered a bit the number of ROIs by using a lax correlation to
%%% see if i get rid of most of the useless things and manage to make this
%%% run better. 


%%% first I need to pick the ones with correlation higher than 0.1

temp_idx_corr=zeros(length(Correlation_group_all_mat),1);
for i=1:length(Correlation_group_all_mat)
    if length(find(Correlation_group_all_mat(:,i)>0.1))>0
     
    temp_idx_corr(i,1)=1;
    else
    temp_idx_corr(i,1)=0;    
    end    
end

idx_corr01=find(temp_idx_corr);

options = statset('UseParallel',1);
[idxKmeans_ZS_CN_corr010_200CL Cmap_ZS_CN_corr010_200CL]=kmeans(ZS_CN(idx_corr01,:),200,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');


%[idxKmeans_ZS_CN_200CL Cmap_ZS_CN_200CL]=kmeans(ZS_CN,200,'Distance','cityblock','Replicates',2,'MaxIter',1000,'Display','final');
%[Model_ZS_CN_200CL,GoodBetas_ZS_CN2_200CL]=Test_Regress(Cmap_ZS_CN,Cmap_ZS_CN,idxKmeans_ZS_CN_200CL,0.1);%%%here we do another linear regression and we select the ones with an r2 value above 0.3


%saveas(gcf,'200rawclusters_CN_f20','jpg');

counter2=0;
for i=1:4
    Fighandle=figure;counter=1;
    set(Fighandle, 'Position', [100, 100, 1200, 1000]);
    for j=1:50
        subplot(7,8,counter);plot(Cmap_ZS_CN_corr010_200CL(j+counter2,:));
        counter=counter+1;
    end
    counter2=counter2+50;
    print(Fighandle,strcat('200clusters_CN_corr010_f20_',num2str(i)),'-dpng','-r0');
    %close all;
end


save ('f20_Kmeans_ZS_CN_corr010_200CL.mat','idxKmeans_ZS_CN_corr010_200CL','Cmap_ZS_CN_corr010_200CL','idx_corr01');


%%

%%% comparing corr and lienar regression results

%%% checking if the ROIs of the rsq are inside the corr



ModelResults_postcorr=[];
parfor i=1:size(ZS_CN(idx_corr),1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregress_CN',ZS_CN(idx_corr(i),:));
    %mdl=stepwiselm(rawregress',ZS(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    %mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_postcorr(i).coef=mdl.Coefficients;
    ModelResults_postcorr(i).MSE=mdl.MSE;
    ModelResults_postcorr(i).Fitted=mdl.Fitted;
    ModelResults_postcorr(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom_corr=[ModelResults_postcorr.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them
idx_corr_rsq=find(rsquare_loom_corr>0.5 & rsquare_loom_corr<1); %%%then select the rsquare that are between 0.3 and 1
figure; %%%to plot them in a raster plot
imagesc(ZS_CN(idx_corr_rsq,:), [-0.5 4]);colormap hot

%saveas(gcf,'rasterplot_r050_f20_CN_3','jpg');

figure; 
plot(mean(ZS_CN(idx_corr_rsq,:)));




%ZS_CN_corr05_rsq05=ZS_CN(idx_corr_rsq,:);
