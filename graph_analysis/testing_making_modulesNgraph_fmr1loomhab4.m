


%%% this script is to analyze the fmr1 loom habituation data with graph theory.
%%% first I will make nodes of the clusters and then make means per genotype.
%%% they I will get the correlation of the responses per loom and make
%%% averages of the fish per genotype. 

%%%% NOTE: the nodes are made base on the
%%% cluster and brain region they belong!!

%%% NOTE2: in this case I am also doing a classification of the ROIs based
%%% on the original wiltype clusters. 

%%%%%%%%%%%%%%%%%%%%%%%%%% i havent finished this code yet. 


%%%from the wildtype loomhab dataset
load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx','S_trim');
 
load('Nodes_N_means_alldatasets.mat','Zbrain_brainMask2D');

%%
load('s20_fmr1_loomhab_CN.mat','ZS_CN');

 load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');

%%% to get the clusters from the Kmeans done to ALL the ROIs
% load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN');

load('s20_good_idx_Fish.mat','idx_Fish');
idx_Fish_cat=categorical(idx_Fish);

load('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','idx_rsq');

load('s20_fmr1_loomhab_CN_part2_High_corr_Nb.mat','High_corr_Nb');

load('fmr1loomhab_BrainRegNclean.mat','PerBrainRegions','RegionList','ROI_temp2','idx_rsq_cleaned');

%% cleaning the clasification of the clusters. 
idx_clean=ismember(idx_rsq,idx_rsq_cleaned);
idx_clean=find(idx_clean);

High_corr_Nb=High_corr_Nb(idx_clean);

%%%% i also need to generate the list of fish saved in the
 %%%% fmr1loomhab_lists.m file
 
 load('s20_fmr1_loomhab_CN_part3.mat','idx_temp1','idx_temp2','idx_temp3','idx_temp4','idx_temp5');
 
 RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

 
 
 %% first clasification based on correlations
 
 
%%% now I will try to see if I can also clasify them with the old clusters
%%% from the wildtype dataset. 

load('rawregressS20.mat') 
load('inhib_s20_regress_200CL.mat','inhib_s20_regress');

rawregress(7,:)=inhib_s20_regress;

figure;
for i=1:size(rawregress,1)
   subplot(2,4,i);
    plot(rawregress(i,:))
end

%%% this one has sound but I only need the other ones. I will also order
%%% them properly from fast hab to inhib. 

rawregressS20(1,:)=rawregress(2,1:904); %%% fasthab sharp
rawregressS20(2,:)=rawregress(4,1:904); %%% fasthab med
rawregressS20(3,:)=rawregress(1,1:904); %%% fasthab broad
rawregressS20(4,:)=rawregress(5,1:904); %%% slopehab
rawregressS20(5,:)=rawregress(3,1:904); %%% nonhab
rawregressS20(6,:)=rawregress(7,1:904); %%% inhib

figure;
for i=1:size(rawregressS20,1)
   subplot(2,4,i);
    plot(rawregressS20(i,:))
end

%%% now i get rid of the first 2 clusters... one is basically in only 2
%%% fish and the other seems to be a drifting effect. 

unique(High_corr_Nb)

figure;
 plot(mean(ZS_CN(idx_rsq_cleaned(find(High_corr_Nb==2)),:)));
 
idx_rsq_good=idx_rsq_cleaned; 

idx_rsq_good(find(High_corr_Nb==1 | High_corr_Nb==2))=[];


%%% to clasify the ROIs with a correlation

Correlation_group_s20={};
counter=1;
for i=1:size(rawregressS20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_CN(idx_rsq_good),1)
        temp_corr=corrcoef(rawregressS20(i,:),ZS_CN(idx_rsq_good(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_s20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_s20=[];
for n=1:size(rawregressS20,1)
Correlation_group_mat_s20(n,:)=Correlation_group_s20{n}(n,:);

end

fmr1_wt_clust_corr=max(Correlation_group_mat_s20);


%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_s20=zeros(length(Correlation_group_mat_s20),1);
for i=1:length(Correlation_group_mat_s20)
    [~,I]=max(Correlation_group_mat_s20(:,i));
    High_corr_Nb_s20(i,1)=I;
    
end

 %%% to check if it worked
 
 figure;
 for i=unique(High_corr_Nb_s20)'
  subplot(2,4,i);   
 plot(mean(ZS_CN(idx_rsq_good(find(High_corr_Nb_s20==i)),:)));
 
 end
 
 
 Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    
    idx_temp=find(High_corr_Nb_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_CN(idx_rsq_good(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_CN(idx_rsq_good(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_rsq_good(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_cat(idx_rsq_good(idx_temp)));%%% for the fish location
    counter=counter+4;
end

 
 %%%% it worked well... although some fasthab clusters dont seem to be very
 %%%% well represented across al fish. 
 
save('fmr1_ROIs_with_WT_clusters.mat','idx_rsq_good','High_corr_Nb_s20','rawregress','fmr1_wt_clust_corr');
 
 %% to get the nodes of the clusters
 
 
  %%% I need to order them as in my draft graphs 
 goodorder_clust=unique(High_corr_Nb_s20);
 
 %%%% merging cluster 3 and 7 (the two broad fasthab) and all the weakly hab (4,5 and 9). 
 

Nodes3=struct;
Nodes3.Mod_loc=[];
Nodes3.Mod_clust=[];
Nodes3.Mod_KmeansID={};
Nodes3.Mod_brain=[];
counter=1;
for clust=goodorder_clust'
       
 idx_temp=idx_rsq_good(find(High_corr_Nb_s20==clust));   
    
 for brain=1:length(RegionList)
     
     
     idx_brain_temp=PerBrainRegions.(RegionList{brain}).idx;
     
     brain_clust_idx=intersect(idx_temp,idx_brain_temp);
 
     if length(brain_clust_idx)<200
       continue 
     elseif length(brain_clust_idx)<500 & length(brain_clust_idx)>200
         moduleN=1;
     elseif length(brain_clust_idx)<1000 & length(brain_clust_idx)>500
     moduleN=2;
     elseif length(brain_clust_idx)<3000 & length(brain_clust_idx)>999
     moduleN=3;
      elseif length(brain_clust_idx)>3000 
     moduleN=4;
     end   
     
options = statset('UseParallel',1); [idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2(brain_clust_idx,:),moduleN,'Options',options,'Distance','sqeuclidean','Replicates',5,'MaxIter',1000,'Display','final');

Nodes3.Mod_loc=vertcat(Nodes3.Mod_loc,Cmap_ROIs);
Nodes3.Mod_clust=vertcat(Nodes3.Mod_clust,ones(size(Cmap_ROIs,1),1)*clust);
Nodes3.Mod_brain=vertcat(Nodes3.Mod_brain,ones(size(Cmap_ROIs,1),1)*brain);


%nodes=find(Nodes2.Mod_clust==clust);
KID=unique(idxKmeans_ROIs);
for ID=1:length(KID)
 Nodes3.Mod_KmeansID{counter,1}=brain_clust_idx(find(idxKmeans_ROIs==ID));   
counter=counter+1;
end
 end
end 

%%% to check if it would work
clust=4;
idx_temp=idx_rsq_good(find(High_corr_Nb_s20==clust)); 
figure;scatter(ROI_temp2(idx_temp,1),ROI_temp2(idx_temp,2));

%%%% i will need to adjust the colors but it seems that it is working. 
figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Nodes.Mod_loc(:,1),ModulNodeses.Mod_loc(:,2),Nodes.Mod_loc(:,3));
gscatter(Nodes3.Mod_loc(:,1),Nodes3.Mod_loc(:,2),Nodes3.Mod_brain);
view(-90,90);

figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Nodes.Mod_loc(:,1),ModulNodeses.Mod_loc(:,2),Nodes.Mod_loc(:,3));
gscatter(Nodes3.Mod_loc(:,1),Nodes3.Mod_loc(:,2),Nodes3.Mod_clust,'ygcbrm','.',20,'on');
view(-90,90);

%%
%%% it seems that it worked... now i need to get ROIs of each fish per
%%% genotype.

fish=vertcat(list1,list2,list3,list4);
list5=union(list1,list3);

%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know
%%% why. might need to check if it is not a mistake in the name. i am
%%% getting rid of it on the list. 
find(idx_Fish==201810048);
fish(find(fish==201810048))=[];


for f=1:length(fish)
    tempfish=find(idx_Fish==fish(f));
    
    if ismember(fish(f),list2)
        group='fmr1';
        %ff=find(ismember(list2,fish(f)));
    elseif ismember(fish(f),list4)
        group='control';
        %ff=find(ismember(list4,fish(f)));
    elseif ismember(fish(f),list5)
        group='hets';
        %ff=find(ismember(list5,fish(f)));
    else 
    end
    temp_group=struct;
  for clust=goodorder_clust'
      
      temp_idx_clust=intersect(idx_rsq_good(find(High_corr_Nb_s20==clust)),tempfish);
     
      
    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes3.Mod_clust==clust);
    for node=(find(Nodes3.Mod_clust==clust))'
    temp_idx=intersect(temp_idx_clust,Nodes3.Mod_KmeansID{node,1});
    
%      figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.
%     figure;
%     plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%     hold on;
%     scatter(ROI_temp2(temp_idx,1),ROI_temp2(temp_idx,2),'filled');
%     view(-90,90);
    
    
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx;  
    
    end

    Nodes3.ROIs_idx.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=idx_ROIs_Node;
    
    
  end
end


 
 %%
 %%% to get the means
 groupnames=fieldnames(Nodes3.ROIs_idx);
 
 
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes3.ROIs_idx.(group));
 for f=1:length(fish)
     
for clust=goodorder_clust'
   

    for node=(find(Nodes3.Mod_clust==clust))'
    
   
    temp_idx=Nodes3.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
%         figure;
%         plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%         hold on;
%         scatter(ROI_temp2(temp_idx,1),ROI_temp2(temp_idx,2));
%         view(-90,90);
    
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        plot(temp_mean);
    end
    
    Nodes3.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

  end
  
  
  %%% making means matrices per fish and testing a correlation
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes3.ROIs_idx.(group));
for f=1:length(fish)

    temp_matrix=NaN(size(Nodes3.Mod_clust,1),904);
    for clust=goodorder_clust'
        
        for node=(find(Nodes3.Mod_clust==clust))'
            
        temp_mean=Nodes3.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes3.mean_matrix.(group).(fish{f})=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes3.corr_matrix.(group).(fish{f})=R_temp;
        
        Nodes3.NaNtest.(group).(fish{f})=isnan(R_temp);
        
end
  end

  %%% testing the matrix with all fish to see if there are still gaps
Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes3.NaNtest.(group).(fish{f}));     
end

Matrix_mean=nanmean(Matrix_mean,3);
  

mean(mean(Matrix_mean))


%% getting the corrmatrix for each loom and each fish
Data_corrMat3=struct;

for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes3.mean_matrix.(group));
    for f=1:length(fish)
    
        temp_mean=Nodes3.mean_matrix.(group).(fish{f});
        
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat3.(group).(fish{f}).loomsR=temp_R;
        
    end
             
end

%% making means of each loom per dataset

for g=1:3
     group=groupnames{g,1};
    
    fish=fieldnames(Data_corrMat3.(group));
    
    Mean_corrMat={};
    for k=1:21
        
        temp_Mean_corrMat=[];
        
        for f=1:length(fish)
            
        temp_mat=Data_corrMat3.(group).(fish{f}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMat3.(group).Mean_corrMat=Mean_corrMat;
    
end


%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes3.Mod_loc(:,1),Nodes3.Mod_loc(:,2),Nodes3.Mod_clust);view(-90,90);
title('Model Nodes');

counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMat3.(group).Mean_corrMat{1,1}); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMat3.(group).Mean_corrMat{1,2}); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMat3.(group).Mean_corrMat{1,3});
subplot(3,8,counter+3);imagesc(Data_corrMat3.(group).Mean_corrMat{1,4});
subplot(3,8,counter+4);imagesc(Data_corrMat3.(group).Mean_corrMat{1,5});
subplot(3,8,counter+5);imagesc(Data_corrMat3.(group).Mean_corrMat{1,6});
subplot(3,8,counter+6);imagesc(Data_corrMat3.(group).Mean_corrMat{1,11});%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMat3.(group).Mean_corrMat{1,12}); %% for 11th loom
title(group);

counter=counter+8;
end

%% now to plot some graphs. 

for g=1:3
    figure;
    count=1;
    group=groupnames{g,1};
   sgtitle(group) 
   
    for k=[2 3 4 5 6 11 12]
    
    R=Data_corrMat3.(group).Mean_corrMat{1,k};
        
        
n=length(Nodes3.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

[~,~,weights] = find(tril(R,-1));

% create the graph object:
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.75; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>

% plot it:
subplot(1,7,count);
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes3.Mod_loc(:,1);
y = Nodes3.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes3.Mod_loc(:,1),Nodes3.Mod_loc(:,2),Nodes3.Mod_clust,'ygcbrm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1

    end
end

%%



%% for the substraction of contros minus fmr1
counter=1;
figure;    
subplot(1,8,counter);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,1},Data_corrMat3.fmr1.Mean_corrMat{1,1})); caxis([-1 1]);colormap('jet');%% for pre loom
subplot(1,8,counter+1);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,2},Data_corrMat3.fmr1.Mean_corrMat{1,2}));caxis([-1 1]); colormap('jet');%% for 1st loom
subplot(1,8,counter+2);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,3},Data_corrMat3.fmr1.Mean_corrMat{1,3}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+3);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,4},Data_corrMat3.fmr1.Mean_corrMat{1,4}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+4);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,5},Data_corrMat3.fmr1.Mean_corrMat{1,5}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+5);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,6},Data_corrMat3.fmr1.Mean_corrMat{1,6}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+6);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,11},Data_corrMat3.fmr1.Mean_corrMat{1,11}));caxis([-1 1]);colormap('jet');%% for 10th loom
subplot(1,8,counter+7);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,12},Data_corrMat3.fmr1.Mean_corrMat{1,12})); caxis([-1 1]);colormap('jet');%% for 11th loom
sgtitle('controls minus fmr1');


%% to check the participation change from 10 to 11th 

%%% this needs to be done with the Brain Connectivity Toolbox. 
%%% so i need to have it in the path. 

W_thr10 = threshold_absolute(abs(Data_corrMat3.fmr1.Mean_corrMat{1,11}), 0.75);

 W_thr11 = threshold_absolute(abs(Data_corrMat3.fmr1.Mean_corrMat{1,12}), 0.75);

[Ppos10 Pneg10] = participation_coef_sign(W_thr10,Nodes3.Mod_clust);

[Ppos11 Pneg11] = participation_coef_sign(W_thr11,Nodes3.Mod_clust);


mean(Ppos10)

ratioPpos10_11=Ppos10./Ppos11;
subsPpos10_11=Ppos10-Ppos11;

figure;histogram(subsPpos10_11);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes3.Mod_loc(:,1),Nodes3.Mod_loc(:,2),25,subsPpos10_11,'filled'); colorbar;colormap('jet')
view(-90,90);
title('fmr1');
%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(Nodes3.Mod_clust(find(subsPpos10_11<-0.3)));

figure;histogram(Nodes3.Mod_clust(find(subsPpos10_11>0.3)));

%% checking a few things


%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for g=1:3
    
  group=groupnames{g,1};
    
    fish=fieldnames(Nodes3.NaNtest.(group));  
     
    
    Matrix_mean=[];
    temp_NaN_nodes=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes3.NaNtest.(group).(fish{f}));  
 temp=double(diag(Nodes3.NaNtest.(group).(fish{f})));
 temp_NaN_nodes=horzcat(temp_NaN_nodes,temp);
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,3,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes3.Mod_loc(:,1),Nodes3.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colorbar;colormap('jet'); caxis([0 1]);
view(-90,90);
title(group);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));

NaN_nodes{g}=temp_NaN_nodes;
end

fish_perNode=[];
meanProp_good=[];
for  g=1:3
    fish_perNode(:,g)=abs(sum(NaN_nodes{g},2)-size(NaN_nodes{g},2));
    meanProp_good(:,g)=1-(sum(NaN_nodes{g},2)/size(NaN_nodes{g},2));
    
end


%%% to see which fish are contributing to specific nodes. 
fish=fieldnames(Nodes3.NaNtest.hets);
fmr1_hindbrain_fish={};
for f=1:length(fish)
    
    temp_mat=Nodes3.NaNtest.hets.(fish{f});
    if temp_mat(26,1)==0 
    fmr1_hindbrain_fish{f}=(fish{f})
    else
    end
end

%%% to see which fish are contributing to specific nodes. 
fish=fieldnames(Nodes3.NaNtest.hets);
fmr1_hindbrain_fish={};
for f=1:length(fish)
    
    temp_mat=Nodes3.NaNtest.hets.(fish{f});
    if temp_mat(26,1)==0 
    fmr1_hindbrain_fish{f}=(fish{f})
    else
    end
end



%% what if i put a threshold on proportion of fish needed to contribute to a node?
%%%% before I run the crosscorrelation. 
%%% with at least 25% of the fish i loose some nodes in the hindbrain although a visual effect still
%%% seems to remain there. 

discard=find(min(meanProp_good,[],2)<0.25); %%% this is the proper way!!! with 0.25 i discard 5 nodes, 0.33=11 nodes 
keep=find(ismember([1:length(Nodes3.Mod_brain)],discard)==0);


%%%% with 0.25, i discarded mostly from contralateral side. 
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes3.Mod_loc(keep,1),Nodes3.Mod_loc(keep,2),Nodes3.Mod_brain(keep));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes3.Mod_loc(:,1),Nodes3.Mod_loc(:,2),Nodes3.Mod_brain(:));view(-90,90);
title('Model Nodes');
a = [1:length(Nodes3.Mod_brain)]'; b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
text(x+dx, y+dy, c);


counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMat3.(group).Mean_corrMat{1,1}(keep,keep)); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMat3.(group).Mean_corrMat{1,2}(keep,keep)); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMat3.(group).Mean_corrMat{1,3}(keep,keep));
subplot(3,8,counter+3);imagesc(Data_corrMat3.(group).Mean_corrMat{1,4}(keep,keep));
subplot(3,8,counter+4);imagesc(Data_corrMat3.(group).Mean_corrMat{1,5}(keep,keep));
subplot(3,8,counter+5);imagesc(Data_corrMat3.(group).Mean_corrMat{1,6}(keep,keep));
subplot(3,8,counter+6);imagesc(Data_corrMat3.(group).Mean_corrMat{1,11}(keep,keep));%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMat3.(group).Mean_corrMat{1,12}(keep,keep)); %% for 11th loom
title(group);

counter=counter+8;
end

%% now to plot some graphs. 

for g=1:3
    figure;
    count=1;
    group=groupnames{g,1};
   sgtitle(group) 
   
    for k=[2 3 4 5 6 11 12]
    
    R=Data_corrMat3.(group).Mean_corrMat{1,k}(keep,keep);
        
        
n=length(Nodes3.Mod_loc(keep));

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

[~,~,weights] = find(tril(R,-1));

% create the graph object:
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.75; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>

% plot it:
subplot(1,7,count);
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes3.Mod_loc(keep,1);
y = Nodes3.Mod_loc(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes3.Mod_loc(keep,1),Nodes3.Mod_loc(keep,2),Nodes3.Mod_clust(keep),'gggbrm','.',20,'off');
hold on;
view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1;

    end
end


%% for the substraction of contros minus fmr1
counter=1;
figure;    
subplot(1,8,counter);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,1}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,1}(keep,keep))); caxis([-1 1]);colormap('jet');%% for pre loom
subplot(1,8,counter+1);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,2}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,2}(keep,keep)));caxis([-1 1]); colormap('jet');%% for 1st loom
subplot(1,8,counter+2);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,3}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,3}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+3);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,4}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,4}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+4);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,5}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,5}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+5);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,6}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,6}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+6);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,11}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,11}(keep,keep)));caxis([-1 1]);colormap('jet');%% for 10th loom
subplot(1,8,counter+7);imagesc(minus(Data_corrMat3.control.Mean_corrMat{1,12}(keep,keep),Data_corrMat3.fmr1.Mean_corrMat{1,12}(keep,keep))); caxis([-1 1]);colormap('jet');%% for 11th loom
sgtitle('controls minus fmr1');


%% to check the participation change from 10 to 11th 
W_thr10 = threshold_absolute(abs(Data_corrMat3.control.Mean_corrMat{1,11}(keep,keep)), 0.75);

 W_thr11 = threshold_absolute(abs(Data_corrMat3.control.Mean_corrMat{1,12}(keep,keep)), 0.75);

[Ppos10 Pneg10] = participation_coef_sign(W_thr10,Nodes3.Mod_clust(keep));

[Ppos11 Pneg11] = participation_coef_sign(W_thr11,Nodes3.Mod_clust(keep));


mean(Ppos10)

ratioPpos10_11=Ppos10./Ppos11;
subsPpos10_11=Ppos10-Ppos11;

figure;histogram(subsPpos10_11);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes3.Mod_loc(keep,1),Nodes3.Mod_loc(keep,2),25,subsPpos10_11,'filled'); colorbar;colormap('jet')
view(-90,90);
title('control');
%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(Nodes3.Mod_clust(find(subsPpos10_11<-0.3)));

figure;histogram(Nodes3.Mod_clust(find(subsPpos10_11>0.3)));


%%
%%%% it seems we will work with the atribution by the kmeans to the nodes. 

save('NodesNgraphFmr1Loomhab3.mat','Nodes3','goodorder_clust','Data_corrMat3','keep','discard');


%%  geting ROIs for unity based on the WT clusters

%%%%


%%% To get the coordinates of each node for unity

for g=1%1:3
     group=groupnames{g,1};
    
     if g==1
         temp_group_idx=idx_temp1%idx_temp5;
     elseif g==2
         temp_group_idx=idx_temp2;
     else
         temp_group_idx=idx_temp4;
     end
     
     
fast_temp_idx1=[];
fast_temp_idx2=[];
for i=1:6

    
    if i==4
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        temp_coor(:,4)=1;
        temp_coor(:,5)=fmr1_wt_clust_corr(temp_idx1(find(ismember(idx_rsq_good(temp_idx1),temp_group_idx))))';
        
        filename=strcat('__Coords_fmr1_clust_slopehab_hab_',group,'_short.csv');        
 
    csvwrite(filename,temp_coor);
    elseif i==5
        temp_idx1=find(High_corr_Nb_s20==i);
        
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        temp_coor(:,4)=1;
        temp_coor(:,5)=fmr1_wt_clust_corr(temp_idx1(find(ismember(idx_rsq_good(temp_idx1),temp_group_idx))))';
        
        filename=strcat('__Coords_fmr1_clust_nonhab_hab_',group,'_short.csv');       
        
    csvwrite(filename,temp_coor);
    elseif i==6
       temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        temp_coor(:,4)=1;
        temp_coor(:,5)=fmr1_wt_clust_corr(temp_idx1(find(ismember(idx_rsq_good(temp_idx1),temp_group_idx))))';
        
        filename=strcat('__Coords_fmr1_clust_inhib_',group,'_short.csv');       
        
    csvwrite(filename,temp_coor);
    
    else       
        temp_idx1=find(High_corr_Nb_s20==i);
        fast_temp_idx1=vertcat(fast_temp_idx1,temp_idx1);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        fast_temp_idx2=vertcat(fast_temp_idx2,temp_idx2);            
    end
           
end
        
        temp_coor=ROI_temp2(fast_temp_idx2,:);
        temp_coor(:,4)=1;
        temp_coor(:,5)=fmr1_wt_clust_corr(fast_temp_idx1(find(ismember(idx_rsq_good(fast_temp_idx1),temp_group_idx))))';
        
filename=strcat('__Coords_fmr1_clust_fasthab_hab_',group,'_short.csv');
csvwrite(filename,temp_coor);

%%% to check the locations
 %figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
 %hold on; scatter(temp_coor(:,1),temp_coor(:,2),'filled');view(-90,90);


end


%%

%%%%%%%%%%%%%%%%  NOTE: there is something weird... I need to check if is
%%%%%%%%%%%%%%%%  an artifact... but it seems that fmr1 fish have far less
%%%%%%%%%%%%%%%%  fasthab ROIs!! 

figure;
counter=1;
for g=[3 1 2]
     group=groupnames{g,1};
    
     if g==1
         temp_group_idx=idx_temp1;%idx_temp5;
     elseif g==2
         temp_group_idx=idx_temp2;
     else
         temp_group_idx=idx_temp4;
     end
         
fast_temp_idx1=[];
fast_temp_idx2=[];
for i=1:6

    subplot(3,6,counter);
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
    hold on;
    
    if i==1
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabSharp_',group,'_ROIs:',num2str(length(temp_coor))));  
        
    elseif i==2
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabMed_',group,'_ROIs:',num2str(length(temp_coor))));    
        
    elseif i==3
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabBroad_',group,'_ROIs:',num2str(length(temp_coor))));    
    
    elseif i==4
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('slopehab_',group,'_ROIs:',num2str(length(temp_coor))));        
    
    elseif i==5
        temp_idx1=find(High_corr_Nb_s20==i);        
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('nonhab_',group,'_ROIs:',num2str(length(temp_coor))));       
           
    elseif i==6
       temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('inhib_',group,'_ROIs:',num2str(length(temp_coor))));              
          
    end
               
    counter=counter+1;
    
end       
        
end



%%% do I have the same number of fasthab ROIs than before?
%%% yes... the difference is only 16.

oldfasthab=(find(High_corr_Nb==3 | High_corr_Nb==7 | High_corr_Nb==8 ));

newfasthab=(find(High_corr_Nb_s20==1 | High_corr_Nb_s20==2 | High_corr_Nb_s20==3 ));

length(oldfasthab)-length(newfasthab)


%%% there is a bigger difference with the slopehabs. 
oldslopehab=(find(High_corr_Nb==6 ));

newslopehab=(find(High_corr_Nb_s20==4 ));

length(oldslopehab)-length(newslopehab)


%%
%%% what if I used the original fmr1 kmeans clusters and sort them as my
%%% main clusters... would I have similar number of ROIs?

figure;
counter=1;
for g=[3 1 2]
     group=groupnames{g,1};
    
     if g==1
         temp_group_idx=idx_temp1;%idx_temp5;
     elseif g==2
         temp_group_idx=idx_temp2;
     else
         temp_group_idx=idx_temp4;
     end
         
fast_temp_idx1=[];
fast_temp_idx2=[];
for i=[8 3 6 4 10]

    subplot(3,5,counter);
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
    hold on;
    
    if i==8 
        temp_idx1=find(High_corr_Nb==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_cleaned(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabSharp_',group,'_ROIs:',num2str(length(temp_coor))));  
        
    elseif i==3  
        temp_idx1=find(High_corr_Nb==i | High_corr_Nb==7);
        temp_idx2=intersect(temp_group_idx,idx_rsq_cleaned(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabMed_',group,'_ROIs:',num2str(length(temp_coor))));    
        
    elseif i==6
        temp_idx1=find(High_corr_Nb==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_cleaned(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('slopehab_',group,'_ROIs:',num2str(length(temp_coor))));        
    
    elseif i==4 
        temp_idx1=find(High_corr_Nb==i | High_corr_Nb==5 | High_corr_Nb==9);        
        temp_idx2=intersect(temp_group_idx,idx_rsq_cleaned(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('nonhab_',group,'_ROIs:',num2str(length(temp_coor))));       
           
    elseif i==10
       temp_idx1=find(High_corr_Nb==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_cleaned(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('inhib_',group,'_ROIs:',num2str(length(temp_coor))));              
          
    end
               
    counter=counter+1;
    
end       
        
end

