
%%% this script is to visualize try to see the distribution of the ROIs each cluster. 
%%% I will be using f20 to compare to the other ones and calculate the
%%% Euclidean distance. 


%%%% distances are in microns according to Gilles


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





%% testing distances

D_fasthab_f60 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,:),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_2_cleaned,:),'euclidean');

M_fasthab = min(D_fasthab_f60');

%%% this is to have a quick look at the descriptive values of the distances
Min_fasthab = min(D_fasthab_f60,[],'all');
Max_fasthab = max(D_fasthab_f60,[],'all');
mean_fasthab = mean(D_fasthab_f60,'all');

%%%% and now actually for the smallest distances
Min_fasthab = min(M_fasthab);
Max_fasthab = max(M_fasthab);
mean_fasthab = mean(M_fasthab);
med_fasthab = median(M_fasthab);
sd_fasthab = std(M_fasthab);
Qs = quantile(M_fasthab,[0.025 0.25 0.50 0.75 0.975]);
figure;histogram(M_fasthab);

% figure;
% scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,2),20,[(M_fasthab./max(M_fasthab))', zeros(numel(M_fasthab),2)],'filled','MarkerFaceAlpha',.2); colorbar;
% map= [(M_fasthab./max(M_fasthab))', zeros(numel(M_fasthab),2)];
% map=sort(map);
% colormap(map);
% 
% figure;
% c=linspace(min(M_fasthab),max(M_fasthab),length(M_fasthab));
% scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,2),20,[(M_fasthab./max(M_fasthab))', zeros(numel(M_fasthab),2)],'filled','MarkerFaceAlpha',.2); colorbar;
% map= [(M_fasthab./max(M_fasthab))', zeros(numel(M_fasthab),2)];
% map=sort(map);
% colormap(map);


% figure;
% c=linspace(min(M_fasthab),max(M_fasthab),length(M_fasthab));
% scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,2),20,[(M_fasthab./max(M_fasthab))', zeros(numel(M_fasthab),2)],'filled','MarkerFaceAlpha',.2); colorbar;
% map= [(c./max(c))', zeros(numel(M_fasthab),2)];
% %map=sort(map);
% colormap(map)


%%% this one looks better!!!
figure;
%c=linspace(min(M),max(M),length(M));
scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,2),20,M_fasthab,'filled','MarkerFaceAlpha',.2); colormap('hot'); colorbar;

%% correlations

%%% checking with its own for f20

ZS_f20=load('f20_CN_r2050_CL5_extra.mat','ZS_CN');

D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,:),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,:),'euclidean');

M_f20=struct;

for i=1:length(D_f20)
        
temp_D=D_f20(i,:);
temp_D(find(temp_D==0))=mean(temp_D);
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0));
M_f20.min(i)=min_temp;
M_f20.idx(i)=idx_temp;

end


for i=1:length(D_f20)
temp_corr=corrcoef(ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned(i),:),ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned(M_f20.idx(i)),:));
Corr_f20(i)=temp_corr(1,2);
end

figure;
scatter(M_f20.min,Corr_f20);

figure;
%c=linspace(min(M),max(M),length(M));
scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,2),20,M_f20.min,'filled','MarkerFaceAlpha',.2); colormap('hot'); colorbar;



%%

%%%% now with f60

ZS_f60=load('final_F60_step1_2.mat','ZS_f60');
load('final_F60_step1_2.mat','ZS_short_F60');


%ZS_f60=load('f60_CN_2.mat','ZS_f60');

%%% to check 
figure;imagesc(ZS_f60.ZS_f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_2_cleaned,ZS_short_F60)); colormap('hot');


%%
D_fasthab_f60 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,:),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_2_cleaned,:),'euclidean');




M_f60 = struct;


for i=1:length(D_fasthab_f60)
        
temp_D=D_fasthab_f60(i,:);
[min_temp idx_temp]=min(temp_D);
M_f60.min(i)=min_temp;
M_f60.idx(i)=idx_temp;

end

for i=1:length(D_fasthab_f60)
temp_corr=corrcoef(ZS_f20.ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned(i),:),ZS_f60.ZS_f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_2_cleaned(M_f60.idx(i)),ZS_short_F60));
Corr_f60(i)=temp_corr(1,2);
end


scatter(M_f60.min,Corr_f60);

figure;
%c=linspace(min(M),max(M),length(M));
scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,2),20,M_fasthab,'filled','MarkerFaceAlpha',.2); colormap('hot'); colorbar;
