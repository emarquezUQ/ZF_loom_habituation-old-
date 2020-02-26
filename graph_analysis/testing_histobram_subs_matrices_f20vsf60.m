
%%%%% this script is to test histograms of the substraction matrices. 
%%%%% I will have a first look at the f20-f60 but the ones I need to check
%%%%% are from the fmr1 dataset. 

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];

cbrewer()

[RdBu]=cbrewer('div','RdBu',101);


%% first the matrices
%%%% to plot the 1st, 10th, 11th and the previous one that look the most
%%%% like 11th
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 data=1;
     datatemp=datasets(data,:);
subplot(2,4,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,4,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,4,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
subplot(2,4,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 4th loom
title(datatemp);
data=2;
     datatemp=datasets(data,:);
subplot(2,4,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,4,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,4,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
subplot(2,4,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,7}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 6th loom
title(datatemp);

%%%%% how doest it look if is not thresholded? 
sub1_weigths=Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep)-Data_corrMat2.(datasets(data+1,:)).Mean_corrMat{1,2}(keep,keep);
sub10_weigths=Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep)-Data_corrMat2.(datasets(data+1,:)).Mean_corrMat{1,11}(keep,keep);
sub11_weigths=Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep)-Data_corrMat2.(datasets(data+1,:)).Mean_corrMat{1,12}(keep,keep);

%%%% is a bit caothic... to many datapoints
figure;
subplot(1,3,1);histogram(sub1_weigths);
subplot(1,3,2);histogram(sub10_weigths);
subplot(1,3,3);histogram(sub11_weigths);

nanmean(sub11_weigths(:))




%%%% this is to get the thresholded data
MatAll_corrected2=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 7 8 9 10 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 6 7 8 9 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected2.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end



counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 data=1;
     datatemp=datasets(data,:);
     loom=fieldnames(MatAll_corrected2.(datatemp));
subplot(2,3,counter);imagesc(MatAll_corrected2.(datatemp).(loom{1}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{1}).Mat); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(MatAll_corrected2.(datatemp).(loom{10}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{10}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(MatAll_corrected2.(datatemp).(loom{11}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{11}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
title(datatemp);


sub1=MatAll_corrected2.(datatemp).(loom{1}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{1}).Mat;
sub10=MatAll_corrected2.(datatemp).(loom{10}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{10}).Mat;
sub11=MatAll_corrected2.(datatemp).(loom{11}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{11}).Mat;


figure;
subplot(1,3,1);histogram(sub1(find(sub1(:))));
subplot(1,3,2);histogram(sub10(find(sub10(:))));
subplot(1,3,3);histogram(sub11(find(sub11(:))));

nanmedian(sub1(find(sub1(:))))
nanmedian(sub10(find(sub10(:))))
nanmedian(sub11(find(sub11(:))))

%%

temp=NaN(length(sub1(:)),2);
temp(1:length(find(sub1(:)>0)),1)=sub1(find(sub1(:)>0));
temp(1:length(find(sub1(:)<0)),2)=abs(sub1(find(sub1(:)<0)));
figure;boxplot(temp);
p = ranksum(x_f20,y_f60);

%%
%%% testing differences without the substraction. i am assuming their distrubution wont be normal
[~, ~, x_f20]=find(MatAll_corrected2.(datatemp).(loom{11}).Mat);
[~, ~, y_f60]=find(MatAll_corrected2.(datasets(data+1,:)).(loom{11}).Mat);

temp=NaN(length(x_f20),2);
temp(1:length(x_f20),1)=x_f20;
temp(1:length(y_f60),2)=y_f60;
figure;boxplot(temp);
p = ranksum(x_f20,y_f60);

%%


