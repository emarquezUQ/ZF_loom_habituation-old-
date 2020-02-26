
%%%% this script is to do the distance and corr analysis. 
%%%% I will be looking at the  distance of between ROIs of each cluster
%%%% and its nearest neighbour. and will also look at how much they
%%%% correlate. I will use as a reference the f20 group. 

%%%% Then, I will look at the distance of the ROIs of each cluster to the
%%%% ROIs of the same cluster in other datasets. Looking for the nearest neighbour. this is
%%%% to see if the distribution is similar (which we expect but we want to quantify). 
%%%% Then I will also look at their correlation. Again, I will use the f20 as reference. 

%%%% NOTE: according to Gilles, the distances are in micros

%%%% thistances are in microns


load('All_More_BrainReg2.mat');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');


datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);


%%% I will start with inside comparison of the f20s


ZS_f20=load('f20_CN_r2050_CL5_extra.mat','ZS_CN');


F20_disNcorr=struct;

for i=1:length(clustersF)

idx_clust=find(strcmp(clustersF, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

%%% to get the distances
D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),'euclidean');

M_f20=struct;

for j=1:length(D_f20)
        
temp_D=D_f20(j,:); %%%% here I am getting each row to be albe to take way the 0 of the ROIs that is being comapared to itself
temp_D(find(temp_D==0))=mean(temp_D); %%% here I am changing the 0 value to the mean of that row. just to be albe to get the min afterwards without getting the 0
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0)); %%% to get the nearest ROI and its location
M_f20.min(j)=min_temp;
M_f20.idx(j)=idx_temp;

end

%%% for the correlation between each ROI and its closest neighbour
Corr_f20=[];
for k=1:length(D_f20)
temp_corr=corrcoef(ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(k),:),ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(M_f20.idx(k)),:));
Corr_f20(k)=temp_corr(1,2);
end


%%% to put all the data together

F20_disNcorr.distance.(char(clustersF(i))).Dis=D_f20;
F20_disNcorr.distance.(char(clustersF(i))).M.min=M_f20.min;
F20_disNcorr.distance.(char(clustersF(i))).M.idx=M_f20.idx;

F20_disNcorr.correlation.(char(clustersF(i))).Corr=Corr_f20;


%%%% to get the descriptive statistics of the distances and corr
F20_disNcorr.distance.(char(clustersF(i))).M.Min_dis = min(M_f20.min);
F20_disNcorr.distance.(char(clustersF(i))).M.Max_dis = max(M_f20.min);
F20_disNcorr.distance.(char(clustersF(i))).M.mean_dis = mean(M_f20.min);
F20_disNcorr.distance.(char(clustersF(i))).M.med_dis = median(M_f20.min);
F20_disNcorr.distance.(char(clustersF(i))).M.sd_dis = std(M_f20.min);
F20_disNcorr.distance.(char(clustersF(i))).M.Qs_dis = quantile(M_f20.min,[0.025 0.25 0.50 0.75 0.975]);


F20_disNcorr.correlation.(char(clustersF(i))).Min_corr = min(Corr_f20);
F20_disNcorr.correlation.(char(clustersF(i))).Max_corr = max(Corr_f20);
F20_disNcorr.correlation.(char(clustersF(i))).mean_corr = mean(Corr_f20);
F20_disNcorr.correlation.(char(clustersF(i))).med_corr = median(Corr_f20);
F20_disNcorr.correlation.(char(clustersF(i))).sd_corr = std(Corr_f20);
F20_disNcorr.correlation.(char(clustersF(i))).Qs_corr = quantile(Corr_f20,[0.025 0.25 0.50 0.75 0.975]);


%clear D_f20 M_f20 Corr_f20


end


%%

%%%% now between f20 and f60. 

ZS_f60=load('final_F60_step1_2.mat','ZS_f60');
load('final_F60_step1_2.mat','ZS_short_F60');

F20vsF60_disNcorr=struct;

for i=1:length(clustersF)

idx_clust=find(strcmp(clustersF, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
temp_fieldname2=fieldnames(f60_cleaned_idxs.clust_f60_CL4_cleaned);temp_fieldname2=temp_fieldname2(idx_clust);


%%% to get the distances
D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname2)),:),'euclidean');

M_f20=struct;

for j=1:size(D_f20,1)
        
temp_D=D_f20(j,:); %%%% here I am getting each row to be albe to take way the 0 of the ROIs that is being comapared to itself
temp_D(find(temp_D==0))=mean(temp_D); %%% here I am changing the 0 value to the mean of that row. just to be albe to get the min afterwards without getting the 0
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0)); %%% to get the nearest ROI and its location
M_f20.min(j)=min_temp;
M_f20.idx(j)=idx_temp;

end

%%% for the correlation between each ROI and its closest neighbour
Corr_f20=[];
for k=1:size(D_f20,1)
temp_corr=corrcoef(ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(k),:),ZS_f60.ZS_f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname2))(M_f20.idx(k)),ZS_short_F60));
Corr_f20(k)=temp_corr(1,2);
end


%%% to put all the data together

F20vsF60_disNcorr.distance.(char(clustersF(i))).Dis=D_f20;
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.min=M_f20.min;
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.idx=M_f20.idx;

F20vsF60_disNcorr.correlation.(char(clustersF(i))).Corr=Corr_f20;


%%%% to get the descriptive statistics of the distances and corr
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.Min_dis = min(M_f20.min);
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.Max_dis = max(M_f20.min);
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.mean_dis = mean(M_f20.min);
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.med_dis = median(M_f20.min);
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.sd_dis = std(M_f20.min);
F20vsF60_disNcorr.distance.(char(clustersF(i))).M.Qs_dis = quantile(M_f20.min,[0.025 0.25 0.50 0.75 0.975]);


F20vsF60_disNcorr.correlation.(char(clustersF(i))).Min_corr = min(Corr_f20);
F20vsF60_disNcorr.correlation.(char(clustersF(i))).Max_corr = max(Corr_f20);
F20vsF60_disNcorr.correlation.(char(clustersF(i))).mean_corr = mean(Corr_f20);
F20vsF60_disNcorr.correlation.(char(clustersF(i))).med_corr = median(Corr_f20);
F20vsF60_disNcorr.correlation.(char(clustersF(i))).sd_corr = std(Corr_f20);
F20vsF60_disNcorr.correlation.(char(clustersF(i))).Qs_corr = quantile(Corr_f20,[0.025 0.25 0.50 0.75 0.975]);


%clear D_f20 M_f20 Corr_f20


end

clear ZS_f60 ZS_short_F60

%%
%%% now for s20

ZS_s20=load('final_S20_step1.mat','ZS_s20');

S_trim=[1:448 453:901 906:1352];

F20vsS20_disNcorr=struct;

for i=1:length(clustersF)

idx_clust=find(strcmp(clustersF, clustersF(i)));    
idx_clust2=find(strcmp(clustersS, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
temp_fieldname2=fieldnames(s20_cleaned_idxs.clust_s20_CL4_cleaned);temp_fieldname2=temp_fieldname2(idx_clust2);


%%% to get the distances
D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname2)),:),'euclidean');

M_f20=struct;

for j=1:size(D_f20,1)
        
temp_D=D_f20(j,:); %%%% here I am getting each row to be albe to take way the 0 of the ROIs that is being comapared to itself
temp_D(find(temp_D==0))=mean(temp_D); %%% here I am changing the 0 value to the mean of that row. just to be albe to get the min afterwards without getting the 0
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0)); %%% to get the nearest ROI and its location
M_f20.min(j)=min_temp;
M_f20.idx(j)=idx_temp;

end

%%% for the correlation between each ROI and its closest neighbour
Corr_f20=[];
for k=1:size(D_f20,1)
temp_corr=corrcoef(ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(k),:),ZS_s20.ZS_s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname2))(M_f20.idx(k)),S_trim));
Corr_f20(k)=temp_corr(1,2);
end


%%% to put all the data together

F20vsS20_disNcorr.distance.(char(clustersF(i))).Dis=D_f20;
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.min=M_f20.min;
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.idx=M_f20.idx;

F20vsS20_disNcorr.correlation.(char(clustersF(i))).Corr=Corr_f20;


%%%% to get the descriptive statistics of the distances and corr
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.Min_dis = min(M_f20.min);
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.Max_dis = max(M_f20.min);
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.mean_dis = mean(M_f20.min);
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.med_dis = median(M_f20.min);
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.sd_dis = std(M_f20.min);
F20vsS20_disNcorr.distance.(char(clustersF(i))).M.Qs_dis = quantile(M_f20.min,[0.025 0.25 0.50 0.75 0.975]);


F20vsS20_disNcorr.correlation.(char(clustersF(i))).Min_corr = min(Corr_f20);
F20vsS20_disNcorr.correlation.(char(clustersF(i))).Max_corr = max(Corr_f20);
F20vsS20_disNcorr.correlation.(char(clustersF(i))).mean_corr = mean(Corr_f20);
F20vsS20_disNcorr.correlation.(char(clustersF(i))).med_corr = median(Corr_f20);
F20vsS20_disNcorr.correlation.(char(clustersF(i))).sd_corr = std(Corr_f20);
F20vsS20_disNcorr.correlation.(char(clustersF(i))).Qs_corr = quantile(Corr_f20,[0.025 0.25 0.50 0.75 0.975]);


%clear D_f20 M_f20 Corr_f20


end

clear ZS_s20 

%%

%%% now for s60

ZS_s60=load('final_S60_step1.mat','ZS_s60');
load('final_S60_step1.mat','ZS_short_S60');

F20vsS60_disNcorr=struct;

for i=1:length(clustersF)

idx_clust=find(strcmp(clustersF, clustersF(i)));    
idx_clust2=find(strcmp(clustersS, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
temp_fieldname2=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname2=temp_fieldname2(idx_clust2);


%%% to get the distances
D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname2)),:),'euclidean');

M_f20=struct;

for j=1:size(D_f20,1)
        
temp_D=D_f20(j,:); %%%% here I am getting each row to be albe to take way the 0 of the ROIs that is being comapared to itself
temp_D(find(temp_D==0))=mean(temp_D); %%% here I am changing the 0 value to the mean of that row. just to be albe to get the min afterwards without getting the 0
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0)); %%% to get the nearest ROI and its location
M_f20.min(j)=min_temp;
M_f20.idx(j)=idx_temp;

end

%%% for the correlation between each ROI and its closest neighbour
Corr_f20=[];
for k=1:size(D_f20,1)
temp_corr=corrcoef(ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(k),:),ZS_s60.ZS_s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname2))(M_f20.idx(k)),ZS_short_S60(S_trim)));
Corr_f20(k)=temp_corr(1,2);
end


%%% to put all the data together

F20vsS60_disNcorr.distance.(char(clustersF(i))).Dis=D_f20;
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.min=M_f20.min;
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.idx=M_f20.idx;

F20vsS60_disNcorr.correlation.(char(clustersF(i))).Corr=Corr_f20;


%%%% to get the descriptive statistics of the distances and corr
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.Min_dis = min(M_f20.min);
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.Max_dis = max(M_f20.min);
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.mean_dis = mean(M_f20.min);
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.med_dis = median(M_f20.min);
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.sd_dis = std(M_f20.min);
F20vsS60_disNcorr.distance.(char(clustersF(i))).M.Qs_dis = quantile(M_f20.min,[0.025 0.25 0.50 0.75 0.975]);


F20vsS60_disNcorr.correlation.(char(clustersF(i))).Min_corr = min(Corr_f20);
F20vsS60_disNcorr.correlation.(char(clustersF(i))).Max_corr = max(Corr_f20);
F20vsS60_disNcorr.correlation.(char(clustersF(i))).mean_corr = mean(Corr_f20);
F20vsS60_disNcorr.correlation.(char(clustersF(i))).med_corr = median(Corr_f20);
F20vsS60_disNcorr.correlation.(char(clustersF(i))).sd_corr = std(Corr_f20);
F20vsS60_disNcorr.correlation.(char(clustersF(i))).Qs_corr = quantile(Corr_f20,[0.025 0.25 0.50 0.75 0.975]);


%clear D_f20 M_f20 Corr_f20


end

clear ZS_s60 ZS_short_S60 S_trim

%% savig

save('distancesNCorr_groups.mat','F20_disNcorr','F20vsF60_disNcorr','F20vsS20_disNcorr','F20vsS60_disNcorr');


%% to graph

%%% to visualize it:


for j=1:length(clustersF)
    
    idx_clust=find(strcmp(clustersF, clustersF(j)));   
    temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

    
    %%%% here i need to change the datasets to visualize
    X=F20_disNcorr.distance.(char(clustersF(j))).M.min';
    Y=F20_disNcorr.correlation.(char(clustersF(j))).Corr';
    numbins = 100;
    
  
    
    [values, centers] = hist3([X Y], [numbins numbins]);
    centers_X = centers{1,1};
    centers_Y = centers{1,2};
    binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
    binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;
    bins_X = zeros(numbins, 2);
    bins_Y = zeros(numbins, 2);
    for i = 1:numbins
        bins_X(i, 1) = centers_X(i) - binsize_X;
        bins_X(i, 2) = centers_X(i) + binsize_X;
        bins_Y(i, 1) = centers_Y(i) - binsize_Y;
        bins_Y(i, 2) = centers_Y(i) + binsize_Y;
    end
    scatter_COL = zeros(length(X), 1);
    onepercent = round(length(X) / 100);
    
    %fprintf('Generating colormap...\n');
    
    for i = 1:length(X)
        if (mod(i,onepercent) == 0)
            fprintf('.');
        end            
        last_lower_X = NaN;
        last_higher_X = NaN;
        id_X = NaN;
        c_X = X(i);
        last_lower_X = find(c_X >= bins_X(:,1));
        if (~isempty(last_lower_X))
            last_lower_X = last_lower_X(end);
        else
            last_higher_X = find(c_X <= bins_X(:,2));
            if (~isempty(last_higher_X))
                last_higher_X = last_higher_X(1);
            end
        end
        if (~isnan(last_lower_X))
            id_X = last_lower_X;
        else
            if (~isnan(last_higher_X))
                id_X = last_higher_X;
            end
        end
        last_lower_Y = NaN;
        last_higher_Y = NaN;
        id_Y = NaN;
        c_Y = Y(i);
        last_lower_Y = find(c_Y >= bins_Y(:,1));
        if (~isempty(last_lower_Y))
            last_lower_Y = last_lower_Y(end);
        else
            last_higher_Y = find(c_Y <= bins_Y(:,2));
            if (~isempty(last_higher_Y))
                last_higher_Y = last_higher_Y(1);
            end
        end
        if (~isnan(last_lower_Y))
            id_Y = last_lower_Y;
        else
            if (~isnan(last_higher_Y))
                id_Y = last_higher_Y;
            end
        end
        scatter_COL(i) = values(id_X, id_Y);
    
    end
    
%     fprintf(' Done!\n');
%     
%     fprintf('Plotting...');
%     
figure;
subplot(1,2,1)
    scatter(X, Y, 10, scatter_COL, 'filled'); colormap('jet'); colorbar;
subplot(1,2,2)  

scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),2),20,X,'filled','MarkerFaceAlpha',.8); colormap('jet'); colorbar;

% subplot(1,3,3)  
% scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),2),20,Y,'filled','MarkerFaceAlpha',.8); colormap('jet'); colorbar;
    
end