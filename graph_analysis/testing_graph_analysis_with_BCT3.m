

%%%%% what if I clean the fish? this script is to redo the analysis but
%%%%% taking into account a minimum number of fish contributing to a node. 



%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for data=1:4
     datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes2.(datatemp).NaNtest_K);  
     
    
    Matrix_mean=[];
    temp_NaN_nodes=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.(datatemp).NaNtest_K.(fish{f}));  
 temp=double(diag(Nodes2.(datatemp).NaNtest_K.(fish{f})));
 temp_NaN_nodes=horzcat(temp_NaN_nodes,temp);
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,4,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colormap('jet'); caxis([0 1]);%colorbar;
view(-90,90);
title(datatemp);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));

NaN_nodes{data}=temp_NaN_nodes;
end

fish_perNode=[];
meanProp_good=[];
for  data=1:4
    fish_perNode(:,data)=abs(sum(NaN_nodes{data},2)-size(NaN_nodes{data},2));
    meanProp_good(:,data)=1-(sum(NaN_nodes{data},2)/size(NaN_nodes{data},2));
    
end

discard=find(min(fish_perNode,[],2)<3); %%% if instead of proportion I use a min of 3 fish? i discard 3 nodes
%discard=find(min(meanProp_good,[],2)<0.25); %%% this is the proper way!!! with 0.25 i discard 1 nodes, 0.33=7 nodes and 0.5=34 nodes
keep=find(ismember([1:102],discard)==0);


%%%% with 0.25, i just discarded 1 of the cerebellum nodes
%%%% with a min of 3 fish i discarded two nodes from the cerebellum and one
%%%% from the tectum. all from the contralateral side. 
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_brain(keep));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain(:));view(-90,90);
title('Model Nodes');
a = [1:102]'; b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
text(x+dx, y+dy, c);


counter=1;
figure;
for data=1:4
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep,keep)); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep));caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep,keep));caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep,keep));caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep));caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep,keep));caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));caxis([-1 1]) %% for 11th loom
title(datatemp);

counter=counter+8;
end


%% testing_graph_analysis_with_BCT3



%%%%% this script is to try to further analyze the graph of the original loomhab data using
%%%%% the Brain Connectivity Toolbox



load('Nodes_N_means_alldatasets2.mat','Nodes2','ROI_temp2_all','Zbrain_brainMask2D');

load('graph_loomNdataset2.mat','Data_corrMat2');

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);


%%% testing the loop works
for data=1:4
     datatemp=datasets(data,:);
        
     fish=fieldnames(Nodes2.(datatemp).corr_matrix);
     
    for tempfish=1:length(fish)
           
    Nodes2.NaNtest.(datatemp).(fish{tempfish})
    end
end

%% checking the corr. distributions
%%%% first just having a look at the corr coeficient values

for data=1:4
    figure;
count=1;
    datatemp=datasets(data,:);
   sgtitle(datatemp) 
   
    for k=[2 3 4 5 6 11 12]
    subplot(1,7,count)
        
    R=Data_corrMat2.(datatemp).Mean_corrMat{1,k}(keep,keep);
        
    [~,~,weights] = find(tril(R,-1));
    %quantile(weights,[0.025 0.25 0.50 0.75 0.975])
    
    count=count+1;
    histogram(weights,20)
    end
    
end



%%
%%%% I will look at density, the number of degrees and strength.
%%%% I could maybe use them to see how they decrease through habituation or
%%%% increase at revovery

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


%kden = density_und(CIJ);

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected.(datatemp).(loom{i}).Mat);

MatAll_corrected.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

%%%% degrees and Strength

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

deg=degrees_und(MatAll_corrected.(datatemp).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected.(datatemp).(loom{i}).Mat));

MatAll_corrected.(datatemp).(loom{i}).deg=deg;
MatAll_corrected.(datatemp).(loom{i}).str=str;
end
end

%%%% to plot histograms
figure;
for i=1:7
subplot(1,7,i);
for data=1:4
datatemp=datasets(data,:);
histogram(MatAll_corrected.(datatemp).(loom{i}).str,20);

hold on;
end
end

%%%% to plot with brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=50; 
    end
    
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).str,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end



%%%% there are many combinations for substractions... depending on what we
%%%% would like to see. i will first test f20 minus f60 and then s20-s60


figure;
counter=1;

for data=3
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    low=-50;high=50;
    
subplot(1,7,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); 
hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,(MatAll_corrected.(datatemp).(loom{i}).deg-MatAll_corrected.(datasets(data+1,:)).(loom{i}).deg),'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end


figure;
counter=1;
for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected.(datatemp).(loom{i}).deg)';
    

end 

subplot(1,4,counter);
imagesc(temp);caxis([0 90])
counter=counter+1;
end


figure;
counter=1;
for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected.(datatemp).(loom{i}).str)';
    

end 

subplot(1,4,counter);
imagesc(temp);caxis([0 60])
counter=counter+1;
end




%% participation and gateway coef

%P = participation_coef(X,modulesID);

%P=participation_coef(Matctrl_1,Nodes2.Mod_clust);
%[Gpos,~]=gateway_coef_sign(Matctrl_1,Nodes2.Mod_clust,1);


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes2.Mod_clust(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected.(datatemp).(loom{i}).Mat,Nodes2.Mod_clust(keep),1);

MatAll_corrected.(datatemp).(loom{i}).P=P;
MatAll_corrected.(datatemp).(loom{i}).Gpos=Gpos;
end
end

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    low=0;high=0.8;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).P,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%%% now with brain regions
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
Pb=participation_coef(temp_mat,Nodes2.Mod_brain(keep));
[Gposb,~]=gateway_coef_sign(MatAll_corrected.(datatemp).(loom{i}).Mat,Nodes2.Mod_brain(keep),1);

MatAll_corrected.(datatemp).(loom{i}).Pb=Pb;
MatAll_corrected.(datatemp).(loom{i}).Gposb=Gposb;
end
end

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    low=0;high=0.8;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).Pb,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end

%% seeing with who they connect

for data=1:4
    datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;

clust_ratio=[];
brain_ratio=[];
clust_ratio_temp=[];
brain_ratio_temp=[];
tot_degrees=[];
clust_deg=[];
brain_deg=[];
clust_deg_temp=[];
brain_deg_temp=[];
for j=1:length(temp_mat)
    temp=find(temp_mat(j,:));
    temp_clust=Nodes2.Mod_clust(keep(temp));
    temp_brain=Nodes2.Mod_brain(keep(temp));
    
    clusters=[4 2 5 6 1 7];%unique(Nodes2.Mod_clust(keep));%%%% in the right order. 
    
    for c=1:length(clusters)
    if isempty(find(temp_clust==clusters(c)))
        clust_ratio_temp(1,c)=0;
        clust_deg_temp(1,c)=0;
    else
        clust_ratio_temp(1,c)=sum(temp_clust==clusters(c))/length(temp); %%% to get the ratio based on all its connections. 
        clust_deg_temp(1,c)=sum(temp_clust==clusters(c));
    end
    end
    
    clust_ratio=vertcat(clust_ratio,clust_ratio_temp);
    clust_deg=vertcat(clust_deg,clust_deg_temp);
    
    
    brains=unique(Nodes2.Mod_brain(keep));
    
    for b=1:length(brains)
    if isempty(find(temp_brain==brains(b)))
        brain_ratio_temp(1,b)=0;
        brain_deg_temp(1,b)=0;
    else
        brain_ratio_temp(1,b)=sum(temp_brain==brains(b))/length(temp); %%% to get the ratio based on all its connections. 
        brain_deg_temp(1,b)=sum(temp_brain==brains(b));
    end
    end
    
    brain_ratio=vertcat(brain_ratio,brain_ratio_temp);
    brain_deg=vertcat(brain_deg,brain_deg_temp);
    
    
    tot_degrees=vertcat(tot_degrees,length(temp));
end


MatAll_corrected.(datatemp).(loom{i}).NodeConnClustRatio=clust_ratio;
MatAll_corrected.(datatemp).(loom{i}).NodeConnBrainRatio=brain_ratio;
MatAll_corrected.(datatemp).(loom{i}).NodeConnDeg=tot_degrees;
MatAll_corrected.(datatemp).(loom{i}).NodeConnClustDeg=clust_deg;
MatAll_corrected.(datatemp).(loom{i}).NodeConnBrainDeg=brain_deg;
end
end

figure;
counter=1;
for data=1:4
    datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=[1 2 3 4 5 6 7]
temp=MatAll_corrected.(datatemp).(loom{i}).NodeConnClustRatio;
subplot(4,7,counter);
imagesc(temp);title(strcat(datatemp,'_',(loom{i})));caxis([0 1]);

counter=counter+1;
    
end
end

figure;
counter=1;
for data=1:4
    datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=[1 2 3 4 5 6 7]
temp=MatAll_corrected.(datatemp).(loom{i}).NodeConnBrainRatio;
subplot(4,7,counter);
imagesc(temp);title(strcat(datatemp,'_',(loom{i})));caxis([0 1]);

counter=counter+1;
    
end
end


%% core vs periphery
% in the function guidelines it says that its supposed to have directed values...
% so i dont know if we can use this, i tried it anyway. 


%CoreP     = core_periphery_dir(W);

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
CoreP=core_periphery_dir(temp_mat);

MatAll_corrected.(datatemp).(loom{i}).CoreP=CoreP;
end
end

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    low=0;high=1;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).CoreP,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%% clustering coef
% i tried it... it might be interesting to see how the connections are lost
% during habituation but I dont think it gives us much more information. 
%%% I could calculate averages and compare between datasets

%Ccoef = clustering_coef_wu(W);
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
Ccoef=clustering_coef_wu(temp_mat);

MatAll_corrected.(datatemp).(loom{i}).Ccoef=Ccoef;
end
end

%%% to test particular combinations. 
data=1;
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));
i=1
h = ranksum(MatAll_corrected.f20.(loom{i}).Ccoef,MatAll_corrected.s20.(loom{i}).Ccoef) %%% 

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    low=0;high=1;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).Ccoef,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%% characteristic path length
%  i first need to get a connection-length matrix, then a distance matrix. I dont think the charpath will
%  really give me something interesting but it can be used to calculate teh
%  small-worldness. I could also compare it between datasets


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
LMat = weight_conversion(temp_mat, 'lengths');
DMat = distance_wei(LMat);
[CharP,~,~,~,~]=charpath(DMat,0,0);

MatAll_corrected.(datatemp).(loom{i}).LMat=LMat;
MatAll_corrected.(datatemp).(loom{i}).DMat=DMat;
MatAll_corrected.(datatemp).(loom{i}).CharP=CharP;
end
end

%%
%%% calculating small-worldness as burgstaller et al 2019 bioRxiv did
%%% although I think I am doing this wrong... I think I should calculate
%%% all the values per fish. 

%%%% so far I can see some differences between genotypes and how the
%%%% differences drop with the stimulus but I am not really sure what to do
%%%% with it. 

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

    SW=mean(MatAll_corrected.(datatemp).(loom{i}).Ccoef)/MatAll_corrected.(datatemp).(loom{i}).CharP;

MatAll_corrected.(datatemp).(loom{i}).SW=SW;

end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).SW;
end
plot(temp);
hold on;

end



%% betweeness centrality. 
%%% i need to first to a connection-length matrix
%%% I tried and the results look kind of weird... not sure what to make of
%%% them. this is supposed to show us some of the main clusters based on
%%% the ones where the shortest paths have to pass through

%BC = betweenness_wei(L);


for data=1:4
datatemp=datasets(data,:);
   
   loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
temp_mat=MatAll_corrected.(datatemp).(loom{i}).LMat;
temp_mat(isnan(temp_mat))=0;    
    
BC=betweenness_wei(temp_mat);
BC=BC/((85-1)*(85-2));%%%% normalizing it
MatAll_corrected.(datatemp).(loom{i}).BC=BC;

end
end

figure;
 counter=1;
for data=1:4
datatemp=datasets(data,:);
   
   
   loom=fieldnames(MatAll_corrected.(datatemp));
for i=[1 2 3 6 7]
 subplot(4,5,counter);   
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).BC,'filled');view(-90,90); colormap('jet');%caxis([0 1]);%colorbar;
title(strcat(datasets(data,:),'_',loom{i}));

counter=counter+1;
end
end



%% modularity and modules based on correlation

% modularity and modules based on correlation
%%% it gives some intresting results but they are not the same as the
%%% HierarchicalConsensus method. it is much fasther though. 


%Ci = modularity_und(W);
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;    
    

[mod, Q]=modularity_und(temp_mat,1);

MatAll_corrected.(datatemp).(loom{i}).mod=mod;
MatAll_corrected.(datatemp).(loom{i}).Q=Q;


end
end

for data=1:4
figure;
%data=1;
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));
for i=1:length(loom)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),MatAll_corrected.(datatemp).(loom{i}).mod);view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
pause(3);
end
end

%%% checking the modularity index
figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).Q;
end
plot(temp);
hold on;

end


figure;
 counter=1;
for data=1:4
datatemp=datasets(data,:);
   
   loom=fieldnames(MatAll_corrected.(datatemp));
for i=[1 2 3 4 5 6 7]
 subplot(4,7,counter);   
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),MatAll_corrected.(datatemp).(loom{i}).mod);view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
title(strcat(datasets(data,:),'_',loom{i}));

counter=counter+1;
end
end

%%

%%%%% HierarchicalConsensus. Based on Betzel 2018 bioRxiv

%%%%%%%%%%%%%%%%    NOTE:
%%%% I havent run this part yet

%%%% I will test with the one that seemed to have worked for the
%%%% fmr1ldataset
%%% I will get some nans to 0 cause I have in some matrices. 
%%% first with f20
temp_mat=MatAll_corrected.f20.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S2 = eventSamples(temp_mat, 150000); % 500000=5 modules. Not sure it was enought though... with 150000 i got the same results. 5 modules
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc2, Tree2] = hierarchicalConsensus(S2,0.05);
[Sall2, thresholds2] = allPartitions(Sc2, Tree2);
%drawHierarchy(Sc2, Tree2)
C2 = coclassificationMatrix(S2);
figure;consensusPlot(C2, Sc2, Tree2)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Sall2(:,1),[],'.',25);view(-90,90);
 

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust(:,1),[],'.',25);view(-90,90);


%%% now with f60
temp_mat=MatAll_corrected.f60.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S3 = eventSamples(temp_mat, 150000); % 150000=10 modules.
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc3, Tree3] = hierarchicalConsensus(S3,0.05);
[Sall3, thresholds3] = allPartitions(Sc3, Tree3);
%drawHierarchy(Sc3, Tree3)
C3 = coclassificationMatrix(S3);
figure;consensusPlot(C3, Sc3, Tree3)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Sall3(:,1),[],'.',25);view(-90,90);

%%% now with s20
temp_mat=MatAll_corrected.s20.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S4 = eventSamples(temp_mat, 150000); % 150000=2 modules. with 500000=2 modules. same ones.
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc4, Tree4] = hierarchicalConsensus(S4,0.05);
[Sall4, thresholds4] = allPartitions(Sc4, Tree4);
%drawHierarchy(Sc4, Tree4)
C4 = coclassificationMatrix(S4);
figure;consensusPlot(C4, Sc4, Tree4)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Sall4(:,1),[],'.',25);view(-90,90);


%%% now with s60
temp_mat=MatAll_corrected.s60.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S5 = eventSamples(temp_mat, 150000); % 150000=4 modules.
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc5, Tree5] = hierarchicalConsensus(S5,0.05);
[Sall5, thresholds5] = allPartitions(Sc5, Tree5);
%drawHierarchy(Sc5, Tree5)
C5 = coclassificationMatrix(S5);
figure;consensusPlot(C5, Sc5, Tree5)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Sall5(:,1),[],'.',25);view(-90,90);



%%% for all of them. 

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
  
    temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
    temp_mat(isnan(temp_mat))=0;

    S = eventSamples(temp_mat, 150000); % 150000 seems to work fine
    %[gamma_min,gamma_max]=gammaRange(temp_mat);
    %S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
    %S = exponentialSamples(temp_mat, 10000);
    [Sc, Tree] = hierarchicalConsensus(S,0.05);
    [Sall, thresholds] = allPartitions(Sc2, Tree);
    %drawHierarchy(Sc2, Tree2)
    C = coclassificationMatrix(S);
    % figure;consensusPlot(C2, Sc2, Tree2)

    % figure;
    % plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
    % hold on;
    % gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Sall2(:,1),[],'.',25);view(-90,90);
    
    
    
MatAll_corrected.(datatemp).(loom{i}).S=S;
MatAll_corrected.(datatemp).(loom{i}).Sc=Sc;
MatAll_corrected.(datatemp).(loom{i}).Tree=Tree;
MatAll_corrected.(datatemp).(loom{i}).Sall=Sall;
MatAll_corrected.(datatemp).(loom{i}).thresholds=thresholds;
MatAll_corrected.(datatemp).(loom{i}).thresholds=thresholds;
MatAll_corrected.(datatemp).(loom{i}).C=C;

end
end


%% substraction of connections
%%%% I am trying f20 minus f60 first. 


f20_f60_subs_mat=struct;   
    
for k=[1:31]
  
    temp_mat1=Data_corrMat2.f20.Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat2=Data_corrMat2.f60.Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat=temp_mat1-temp_mat2;
    
    f20_f60_subs_mat.Subs_Mean_corrMat{1,k}=temp_mat;
    
end


f20_s20_subs_mat=struct;   
    
for k=[1:31]
  
    temp_mat1=Data_corrMat2.f20.Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat2=Data_corrMat2.s20.Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat=temp_mat1-temp_mat2;
    
    f20_s20_subs_mat.Subs_Mean_corrMat{1,k}=temp_mat;
    
end
%%%%% and doing a ratio based on first loom response first
%%%% there a few connections that increase but I am not sure if they are
%%%% meaningful 

for data=1:4
datatemp=datasets(data,:);
for k=[1:31]
  
    temp_mat1=Data_corrMat2.(datatemp).Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat2=Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep); %%% to get the first loom
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat=temp_mat1./temp_mat2;
    
    f20_f60_subs_mat.(datatemp).RatiosLoom1{1,k}=temp_mat;
    
    %%% to clean based on the previous threshold
    temp_mat1_idx=find(threshold_absolute(abs(temp_mat1),0.75));
    %temp_mat2_idx=find(threshold_absolute(abs(temp_mat2),0.75));
    %temp_mat_idx=union(temp_mat1_idx,temp_mat2_idx);
    
    temp_mat_cleaned=zeros(size(temp_mat));
    temp_mat_cleaned(temp_mat1_idx)=temp_mat(temp_mat1_idx);
    
    f20_f60_subs_mat.(datatemp).RatiosLoom1_cleaned{1,k}=temp_mat_cleaned;
end
end
%%
%%%% here i am trying to get the connections from my previous filtering
  
for k=[1:31]
  
    temp_mat1=Data_corrMat2.f20.Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat1_idx=find(threshold_absolute(abs(temp_mat1),0.75));
    
    temp_mat2=Data_corrMat2.f60.Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat2_idx=find(threshold_absolute(abs(temp_mat2),0.75));
    
    temp_mat=temp_mat1-temp_mat2;
    
    temp_mat_idx=union(temp_mat1_idx,temp_mat2_idx);
    
    temp_mat_cleaned=zeros(size(temp_mat));
    temp_mat_cleaned(temp_mat_idx)=temp_mat(temp_mat_idx);
    
    f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,k}=temp_mat_cleaned;
    
end


for k=[1:31]
  
    temp_mat1=Data_corrMat2.f20.Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat1_idx=find(threshold_absolute(abs(temp_mat1),0.75));
    
    temp_mat2=Data_corrMat2.s20.Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat2_idx=find(threshold_absolute(abs(temp_mat2),0.75));
    
    temp_mat=temp_mat1-temp_mat2;
    
    temp_mat_idx=union(temp_mat1_idx,temp_mat2_idx);
    
    temp_mat_cleaned=zeros(size(temp_mat));
    temp_mat_cleaned(temp_mat_idx)=temp_mat(temp_mat_idx);
    
    f20_s20_subs_mat.Subs_Mean_corrMat_cleaned{1,k}=temp_mat_cleaned;
    
end
    
    %% now to plot some graphs. 
    
%     figure;
%      count=1;
%     for data=1:4
%  datatemp=datasets(data,:);
 figure;
%sgtitle(datatemp);
 count=1;
    for k=[2 3 4 5 6 11 12]
   % figure;
    %set(gcf, 'Position',  [200, 200, 700, 900]);
    set(gcf, 'Position',  [200, 200, 1200, 900]);
    
    %R=Data_corrMat2.(datatemp).Mean_corrMat{1,k};
    %R=f20_f60_subs_mat.Subs_Mean_corrMat{1,k};
    %R=f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,k}; 
    R=f20_s20_subs_mat.Subs_Mean_corrMat_cleaned{1,k};
    %R=f20_f60_subs_mat.(datatemp).RatiosLoom1_cleaned{1,k}; 
    
n=length(Nodes2.Mod_loc(keep));

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%[~,~,weights] = find(tril(R,-1));
weights = nonzeros(tril(R,-1));

%quantile(abs(weights),[0.025 0.25 0.50 0.75 0.975])


% create the graph object:
%G = graph(s,t,weights,n);
G = graph(R);

% mark the lines to remove from the graph:
threshold = 0; %  minimum correlation to plot
%threshold = 0; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>


subplot(1,7,count);
% plot it:
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
  colormap jet;caxis([-1 1]);%colorbar;
  %colormap jet;caxis([0 3]);%colorbar;
  p.EdgeCData=G.Edges.Weight;
%p.EdgeColor = [G.Edges.Weight>0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight<0.']; %% red high, blue low
%p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];

% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*3;
%p.LineWidth = 1;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes2.Mod_loc(keep,1);
y = Nodes2.Mod_loc(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
%gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust,'rgggbm','.',20,'off');
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust(keep),'ygcbrm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

%saveas(gcf,strcat('subsf20vsf60_',num2str(k),'.svg'));

count=count+1;

    end

    %end

   