%%% for f20

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

load('f20_r2050_CL5.mat')


goodmaps=[11 31 38 39 41 42]; %%% i took the sound response out
for i=1:length(goodmaps)
rawregress(i,:)=Cmap_ZS(goodmaps(i),:);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(goodmaps);

for i=1:length(goodmaps)
   
    subplot(3,2,counter);plot(rawregress(i,:));
     counter=counter+1;
end

print(Fighandle,'goodmaps_F20_poster','-dpdf','-bestfit');


%%%this is to make a figure where we will plot the mean of clusters
%%%selected the raster plot, the levels and the fish where they are found.
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(goodmaps);counter=1;
for i=goodmaps
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans_ZS==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_temp,:),1)); ylim([-2 inf])%%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_temp,:),[0 3]); %%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_temp));%%% for the fish location
    counter=counter+4;
end

print(Fighandle,'raw_multigraph_F20_poster2','-dpdf','-bestfit');


%%% for the raster plots
imagesc(ZS, [-0.5 4]);colormap hot

print(gcf,'rasterplot_allZS_f20_poster','-dpdf','-bestfit');

imagesc(ZS(idx_rsq,:), [-0.5 4]);colormap hot

print(gcf,'rasterplot_r050_f20_poster','-dpdf','-bestfit');


GoodBetas_ZS2=[1 3 4];


idxKmeans_ZS_rsq_poster=idxKmeans_ZS_rsq;
idxKmeans_ZS_rsq_poster(find(idxKmeans_ZS_rsq==5))=4;

idxKmeans1_coef_rsq=idxKmeans_ZS_rsq_poster;%%% to make a variable with the common idx of the interesting clusters from the 1st kmeans and the traces hat passed the regression
rows=length(GoodBetas_ZS2);
counter=1;


idxKmeans_final=zeros(size(ZS,1),1);%%%to make a variable with 0s the size of ZS
idxKmeans_final(idx_coef_rsq)=idxKmeans1_coef_rsq;%%% and put the filtered responses by r2 and filters

temp=[];
counter=1;
for i=GoodBetas_ZS2 %%%to take the clusters we want
    idx_temp=find(idxKmeans_final==i);   %%%and put the filtered values in a new variable in different columns
    temp{counter}=idx_temp;    
    counter=counter+1;    
end

Numbers(1)=0; %%% to change the first value of Numbers to 0 (cause it was not relevant)
%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(GoodBetas_ZS2),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
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
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_ZS2)));yplot=ceil(length(GoodBetas_ZS2)/xplot);
for i=GoodBetas_ZS2
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1),'color',colors(counter2,:)/256); ylim([-2 inf]);%%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end



print(gcf,'good_multigraph_r050_CL3_F20_poster2','-dpdf','-bestfit');



%%% for a sample fish, i will try with fish 13

idx_fish13_filtered=find(idx_Fish(idx_coef_rsq)==13);


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_ZS2)));yplot=ceil(length(GoodBetas_ZS2)/xplot);
for i=GoodBetas_ZS2
    
    
    idx_temp=intersect(find(idxKmeans1_coef_rsq==i),idx_fish13_filtered);
    
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end



print(Fighandle,'good_multigraph_r050_CL3_fish13_F20_poster','-dpdf','-bestfit');


save('f20_r2050_CL3_idxs_poster.mat','MatFiles','idx_Fish','idx_Plane','idx_coef_rsq','idxKmeans1_coef_rsq','GoodBetas_ZS2','-v7.3');


%%

%%% for f60


load('f60_r2050_CL4.mat')

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(goodmaps);counter=1;
for i=goodmaps
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1));ylim([-2 inf]); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end

print(Fighandle,'raw_multigraphs_f60_poster2','-dpdf','-bestfit');

%%

%%% for s20
load('s20_r2050_CL6.mat')

goodmaps=[12    19    26    30    41];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(goodmaps);counter=1;
for i=goodmaps
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1));ylim([-2 inf]); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end

print(Fighandle,'raw_multigraphs_S20_poster2','-dpdf','-bestfit');


%%

%%% for s60


load('s60_r2050_CL3.mat')

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(goodmaps);counter=1;
for i=goodmaps
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1));ylim([-2 inf]); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end

print(Fighandle,'raw_multigraphs_S60_poster2','-dpdf','-bestfit');

save('s60_r2050_CL3_idxs_poster.mat','MatFiles','idx_Fish','idx_Plane','idx_coef_rsq','idxKmeans1_coef_rsq','goodmaps','-v7.3');

%%% for the coloured multigraph

goodmaps=[8 37 3];

%%%%%this is to do it together with the raster plots and histograms of
%%%%%localization of the clusters in fish and levels
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(goodmaps)));yplot=ceil(length(goodmaps)/xplot);
for i=goodmaps
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1),'color',colors(counter2,:)/256);ylim([-2 inf]);xlim([-200 3700]); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
     subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
     subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    %counter=counter+2;
    counter2=counter2+1
     counter=counter+4;
end

print(Fighandle,'good_multigraphs_S60_poster3','-dpdf','-bestfit');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(goodmaps)));yplot=ceil(length(goodmaps)/xplot);
for i=goodmaps
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,2,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1),'color',colors(counter2,:)/256,'LineWidth',2);ylim([-2 inf]);xlim([-200 3700]); %%%to plot the mean
    subplot(rows,2,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);xlim([-200 3700]);%%%for the raster plot
%      subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
%      subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+2;
    counter2=counter2+1
     %counter=counter+4;
end

print(Fighandle,'good_multigraphs_S60_poster4_2','-dpdf','-bestfit');


